import openmm as mm
import openmm.app as app
import openmm.unit as unit
import pandas as pd
import numpy as np
import calvados
import os
import argparse
import mpipi_GG as mpipi 
import pandas as pd

def distribute_points_hexagonal(N, l1):
    """Distributes N points evenly over a square plane using a hexagonal lattice."""
    aspect_ratio = np.sqrt(3) / 2  # Height to width ratio in hex grid
    num_cols = int(np.ceil(np.sqrt(N / aspect_ratio)))
    num_rows = int(np.ceil(N / num_cols))
    
    separation_x = l1 / (num_cols - 0.5)
    separation_y = l1 / (num_rows * aspect_ratio)
    radius = min(separation_x / 2, separation_y / (2 * np.sqrt(3) / 3))
    
    points = []
    for row in range(num_rows):
        y = row * separation_y
        x_offset = separation_x / 2 if row % 2 else 0
        num_cols_in_row = num_cols - 1 if row % 2 else num_cols
        for col in range(num_cols_in_row):
            x = x_offset + col * separation_x
            if x <= l1 and y <= l1:
                points.append((x, y))
                if len(points) == N:
                    break
        if len(points) == N:
            break
    return points, radius

def generate_data(proteins, Nchains, sigma, starting_separation_between_residues, dataset,lz):
    """Generates the data array with amino acids, positions, and chain information."""
    # Get the sequence and side length
    sequence = dataset['Sequence'][proteins]
    l1 = len(sequence)
    
    # Generate starting positions
    start_positions, _ = distribute_points_hexagonal(Nchains, l1)
    
    # Initialize variables
    data_list = []
    atom_counter = 0  # Atom counter
    
    # Build atoms
    for chain_idx, (x_start, y_start) in enumerate(start_positions):
        chain_id = chain_idx + 1
        zstart = 0.5*lz+np.random.normal(0.0, sigma, 1)[0]  # Starting z-position for this chain
        
        # Generate atoms for this chain
        for i, amino_acid in enumerate(sequence):
            atom_counter += 1
            # Map amino acid to index in HPSUrry.amino_acids
            amino_acid_index = calvados.reverse[calvados.reverse_single[amino_acid]]  # This gives the index
            x = x_start
            y = y_start
            z = zstart + (i * starting_separation_between_residues -
                          (0.5 * starting_separation_between_residues * len(sequence)))
            
            # Append atom to the list
            data_list.append([atom_counter, chain_id, amino_acid_index, x, y, z])
    
    # Convert data_list to numpy array
    data = np.array(data_list)
    return data

def create_topology(data):
    """Creates the OpenMM topology from the data array."""
    topology = app.Topology()
    current_chain_id = None
    atom1 = None
    chain = None

    for i in range(len(data)):
        chain_id = data[i, 1]
        amino_acid_index = data[i, 2]
        amino_acid_name = calvados.amino_acids[amino_acid_index]


        if chain_id != current_chain_id:
            current_chain_id = chain_id
            chain = topology.addChain()
            residue = topology.addResidue(name="CG-residue", chain=chain)
            atom1 = topology.addAtom(name=amino_acid_name, element=None, residue=residue)
        else:
            residue = topology.addResidue(name="CG-residue", chain=chain)
            atom2 = topology.addAtom(name=amino_acid_name, element=None, residue=residue)
            topology.addBond(atom1, atom2)
            atom1 = atom2
    return topology

def create_system_calvados(topology, data, pbc, lxy, lz, temperature):
    """Creates the OpenMM system using the Calvados model."""
    # Convert positions to nanometers
    positions = (data[:, 3:6] / 10.0) * unit.nanometer

    # Initialize system
    system = mm.System()

    # Add particles with mass adjustments
    for chain in topology.chains():
        num_a = chain.topology.getNumAtoms()
        counted_a = 0
        for atom in chain.atoms():
            if counted_a == 0:
                factor = 2
            elif counted_a == num_a-1:
                factor = 16
            else:
                factor=0
            system.addParticle(calvados.mass[atom.name]+factor)
            counted_a+=1


    ### MAKE BONDS ###
    hbond = mm.HarmonicBondForce()
    for bond in topology.bonds():
        hbond.addBond(bond.atom1.index, bond.atom2.index, 0.38, 2000*calvados._kcal_to_kj)
    
    hbond.setForceGroup(1)

    energy_expression = 'select(step(r-2^(1/6)*s),4*eps*l*((s/r)^12-(s/r)^6-shift),4*eps*((s/r)^12-(s/r)^6-l*shift)+eps*(1-l))'
    ah = mm.CustomNonbondedForce(energy_expression+'; s=0.5*(s1+s2); l=0.5*(l1+l2); shift=(0.5*(s1+s2)/rc)^12-(0.5*(s1+s2)/rc)^6')
    
    ah.addGlobalParameter('eps',calvados.epsilon*unit.kilojoules_per_mole)
    cutoff = 2.0
    ah.addGlobalParameter('rc',cutoff*unit.nanometer)
    ah.addPerParticleParameter('s')
    ah.addPerParticleParameter('l')
    
    print('rc',cutoff*unit.nanometer)
    kT = 8.3145*temperature*1e-3
    fepsw = lambda T : 5321/T+233.76-0.9297*T+0.1417*1e-2*T*T-0.8292*1e-6*T**3
    epsw = fepsw(temperature)
    lB = 1.6021766**2/(4*np.pi*8.854188*epsw)*6.022*1000/kT
    ionic = 0.150
    yukawa_kappa = np.sqrt(8*np.pi*lB*ionic*6.022/10)
    yu = mm.CustomNonbondedForce('q*(exp(-kappa*r)/r-shift); q=q1*q2')
    yu.addGlobalParameter('kappa',yukawa_kappa/unit.nanometer)
    yu.addGlobalParameter('shift',np.exp(-yukawa_kappa*4.0)/4.0/unit.nanometer)
    yu.addPerParticleParameter('q')
    
    for chain in topology.chains():
        num_a = chain.topology.getNumAtoms()
        counted_a = 0
        for atom in chain.atoms():
            if counted_a == 0:
                factor = 1
            elif counted_a == num_a-1:
                factor = -1
            else:
                factor=0
            yu.addParticle([(calvados.charge[atom.name]+factor)*np.sqrt(lB*kT)])
            ah.addParticle([calvados.size[atom.name]/10*unit.nanometer, calvados.hydropathy[atom.name]*unit.dimensionless])
            counted_a +=1
    yu.createExclusionsFromBonds([(bond[0].index, bond[1].index) for bond in topology.bonds()], 1)
    ah.createExclusionsFromBonds([(bond[0].index, bond[1].index) for bond in topology.bonds()], 1)
    yu.setForceGroup(2)
    ah.setForceGroup(3)
    yu.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
    ah.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
    hbond.setUsesPeriodicBoundaryConditions(True)
    yu.setCutoffDistance(4*unit.nanometer)
    ah.setCutoffDistance(cutoff*unit.nanometer)
    
    ### REMOVE CMM MOTION ###
    force = mm.CMMotionRemover()
    force.setForceGroup(5)

    ### GENERATE SYSTEM ###
    system.addForce(hbond)
    system.addForce(ah)
    system.addForce(yu)
    system.addForce(force)

    box_vec_a = np.array([lxy, 0, 0])*unit.nanometer
    box_vec_b = np.array([0, lxy, 0])*unit.nanometer
    box_vec_c = np.array([0, 0, lz])*unit.nanometer
    system.setDefaultPeriodicBoxVectors(box_vec_a, box_vec_b, box_vec_c)

    return system, positions

def run_simulation(system, topology, positions, experiment_name, runtime, temperature, sampling):
    """Runs the OpenMM simulation with specified parameters."""
    # Write initial PDB file
    with open(f'{experiment_name}.pdb', 'w') as f:
        app.PDBFile.writeFile(topology, positions, f)
    
    # Set up integrator
    temperature_unit = temperature * unit.kelvin
    friction_coeff = 1 / unit.picosecond
    timestep = 10 * unit.femtosecond
    integrator = mm.LangevinMiddleIntegrator(temperature_unit, friction_coeff, timestep)
    
    # Set up platform
    properties = {'Precision': 'mixed'}
    platform_name = 'CUDA'
    platform = mm.Platform.getPlatformByName(platform_name)
    
    # Set up simulation
    simulation = app.Simulation(topology, system, integrator, platform, properties)
    simulation.context.setPositions(positions)
    simulation.minimizeEnergy()
    simulation.context.setVelocitiesToTemperature(temperature_unit)
    
    # Set up reporters
    sampling_interval = int(sampling)
    simulation.reporters.append(app.DCDReporter(f'{experiment_name}.dcd', sampling_interval, enforcePeriodicBox=True))
    simulation.reporters.append(app.StateDataReporter(f'{experiment_name}.csv', sampling_interval, step=True, time=True,
                                                      potentialEnergy=True, kineticEnergy=True, totalEnergy=True,
                                                      temperature=True, speed=True))
    
    # Run simulation
    if runtime == 0:
        simulation.runForClockTime(20 * unit.hours)
    else:
        simulation.step(int(runtime))
    
    # Save checkpoint and state
    simulation.saveCheckpoint(f'{experiment_name}.chk')
    state = simulation.context.getState(getPositions=True, getVelocities=True, getForces=True,
                                        enforcePeriodicBox=True)
    with open(f'{experiment_name}_state.xml', 'w') as f:
        f.write(mm.XmlSerializer.serialize(state))


def rerun(system, topology, positions, experiment_name, runtime, temperature, sampling):
        with open(f'{experiment_name}.pdb', 'w') as f:
            app.PDBFile.writeFile(topology, positions, f)
    
        # Set up integrator
        temperature_unit = temperature * unit.kelvin
        friction_coeff = 1 / unit.picosecond
        timestep = 10 * unit.femtosecond
        integrator = mm.LangevinMiddleIntegrator(temperature_unit, friction_coeff, timestep)
        
        # Set up platform
        properties = {'Precision': 'mixed'}
        platform_name = 'CUDA'
        platform = mm.Platform.getPlatformByName(platform_name)
        
        # Set up simulation
        simulation = app.Simulation(topology, system, integrator, platform, properties)
        simulation.context.setPositions(positions)
        simulation.minimizeEnergy()
        simulation.context.setVelocitiesToTemperature(temperature_unit)

        ### REPORTERS ###
        dcd_reporter = app.DCDReporter(f'{experiment_name}.dcd', sampling, enforcePeriodicBox=True,append=True)
        state_data_reporter = app.StateDataReporter(f'{experiment_name}.csv', sampling, step=True, time=True, potentialEnergy=True,
                                                    kineticEnergy=True, totalEnergy=True, temperature=True, speed=True,append=True)
        simulation.reporters.append(dcd_reporter)
        simulation.reporters.append(state_data_reporter)
        
        ### RUN SIMULATION ###
        
        simulation.loadCheckpoint(f'{experiment_name}.check')
        if runtime == 0:
            simulation.runForClockTime(20)
        else:
            simulation.step(runtime)
        simulation.saveCheckpoint(f'{experiment_name}.check')
        state = simulation.context.getState(getPositions=True, getVelocities=True, getForces=True,enforcePeriodicBox=True)
        with open(f'{experiment_name}_prev.xml', 'w') as f:
            f.write(mm.XmlSerializer.serialize(state))

def run_calvados(temperature,protein_num,filename,Nchains,sigma,lxy,lz,sampling):
    # Parameters
    starting_separation_between_residues = 3.2
    pbc = True
    runtime = 0  # Set to 0 to run for clock time

    # Read the dataset
    dataset = pd.read_csv(filename)
    experiment_name = dataset['Protein Name'][protein_num]
    data = generate_data(protein_num, Nchains, sigma, starting_separation_between_residues, dataset, lz)
    topology = create_topology(data)
    system, positions = create_system_calvados(topology, data, pbc, lxy, lz, temperature)

    if os.path.isfile(f'{experiment_name}.check'):
        rerun(system, topology, positions, experiment_name, runtime, temperature, sampling)
    else:
        run_simulation(system, topology, positions, experiment_name, runtime, temperature, sampling)

def generate_lammps_input(sampling, temperature, protein_name):
    """Generates a LAMMPS input script with the given sampling, temperature, and data file name."""
    with open("lammps.in", "w") as file:
        file.write(f"""variable dt equal 10
variable temp equal {temperature}
variable sample equal {int(sampling)}

units real
dimension 3
boundary p p p
processors * * *
atom_style full
bond_style harmonic
dielectric 80.0

read_data {protein_name}.lmp
include pairs.txt
special_bonds fene

neighbor 3.5 multi
neigh_modify every 10 delay 0
comm_style tiled


## ########################################################
##      PART 2: Prepare for the actual simulation
## ########################################################

timestep ${{dt}}
timer normal

fix fxnve all nve
fix fxlange all langevin ${{temp}} ${{temp}} 5000.0 1234
fix fxbal all balance 1000 1.05 rcb


## ########################################################
##      PART 3: SETUP OUTPUT INFO
## ########################################################

dump 1 all custom ${{sample}} traj.dump mol type x y z
dump_modify 1 sort id

thermo ${{sample}}
thermo_style custom step pe press ke lz temp spcpu
thermo_modify flush yes


## ########################################################
##      PART 4: INITIALIZE SYSTEM
## ########################################################

velocity all create ${{temp}} 1234 rot yes dist gaussian
minimize 0.0 0.00000001 1000 100000
timer timeout 22:30:00
run 1000000000
write_restart restart.equil
""")

def generate_lammps_restart_input(sampling, temperature):
    """Generates a LAMMPS input script to restart from restart.equil if available."""

    # Generate the restart input script
    with open("lammps_restart.in", "w") as file:
        file.write(f"""variable dt equal 10
variable temp equal {temperature}
variable sample equal {int(sampling)}

units real
dimension 3
boundary p p p
processors * * *
atom_style full
bond_style harmonic
dielectric 80.0

# Read from restart file instead of data file
read_restart restart.equil
include pairs.txt
special_bonds fene

neighbor 3.5 multi
neigh_modify every 10 delay 0
comm_style tiled


## ########################################################
##      PART 2: Prepare for the actual simulation
## ########################################################

timestep ${{dt}}
timer normal

fix fxnve all nve
fix fxlange all langevin ${{temp}} ${{temp}} 5000.0 1234
fix fxbal all balance 1000 1.05 rcb


## ########################################################
##      PART 3: SETUP OUTPUT INFO
## ########################################################

dump 1 all custom ${{sample}} traj_restart.dump mol type x y z
dump_modify 1 sort id

thermo ${{sample}}
thermo_style custom step pe press ke lz temp spcpu
thermo_modify flush yes


## ########################################################
##      PART 4: RESUME SYSTEM FROM RESTART
## ########################################################

velocity all create ${{temp}} 1234 rot yes dist gaussian
timer timeout 22:30:00
run 1000000000
write_restart restart.equil
""")


def run_mpipi(temperature,protein_num,filename,Nchains,sigma,lxy,lz,sampling):
    lxy *= 10
    lz *= 10
    dataset = pd.read_csv(filename)
    sequence = dataset['Sequence'][protein_num]
    experiment_name = dataset['Protein Name'][protein_num]
    starting_separation_between_residues = 0.7
    start_positions, radius = distribute_points_hexagonal(Nchains, lxy)
    atoms_list = []
    bonds = []
    atom_counter = 1  # Global atom counter
    bond_counter = 1  # Global bond counter

    # Build atoms and bonds
    for chain_idx, (x_start, y_start) in enumerate(start_positions):
        zstart = np.random.normal(0.0, sigma, 1)[0]  # Starting z-position for this chain
        chain_atoms = []

        # Generate atoms for this chain
        for i, amino_acid in enumerate(sequence):
            # Atom properties
            atom_number = atom_counter
            chain_id = chain_idx + 1
            amino_acid_type = mpipi.amino_acids[mpipi.reverse_single[amino_acid]]
            x = x_start
            y = y_start
            z = zstart + (i * starting_separation_between_residues - (0.5 * starting_separation_between_residues * len(sequence)))

            # Append atom to the list
            atoms_list.append([atom_number, chain_id, amino_acid_type, 0, x, y, z])
            chain_atoms.append(atom_number)
            atom_counter += 1

        # Generate bonds for this chain
        for i in range(len(chain_atoms) - 1):
            bonds.append([bond_counter, 1, chain_atoms[i], chain_atoms[i + 1]])  # Assuming bond type 1
            bond_counter += 1

    # Convert atoms list to DataFrame
    atoms = pd.DataFrame(atoms_list, columns=['Atom Number', 'Chain', 'Amino Acid', 'Ignore', 'x', 'y', 'z'])
    if not os.path.isfile("restart.equil"):
        generate_lammps_input(sampling, temperature, experiment_name)
    else:
        generate_lammps_restart_input(sampling, temperature)
    # Now write the data to the LAMMPS file
    with open(f'{experiment_name}.lmp', 'w') as file:
        # Write header
        file.write("LAMMPS Description file\n\n")
        file.write(f"{len(atoms)} atoms\n")
        file.write(f"{len(bonds)} bonds\n")
        file.write(f"{len(mpipi.amino_acids) * 2 + 4} atom types\n")
        file.write(f"2 bond types\n\n")
        
        # Box dimensions
        file.write(f"{0} {lxy} xlo xhi\n")
        file.write(f"{0} {lxy} ylo yhi\n")
        file.write(f"{0} {lz} zlo zhi\n\n")
        
        rna_mass = [329.2, 305.2, 345.2, 306.2]
        # Masses
        file.write("Masses\n\n")
        scale = [1, 0.7]
        for i in range(2):
            for aa in range(len(mpipi.amino_acids)):
                mass = scale[i] * mpipi.mass[mpipi.reverse[aa + 1]]
                file.write(f"{(aa + 1) + (i * 20)} {mass:.3f}\n")
        for rm in range(len(rna_mass)):
            file.write(f"{41 + rm} {rna_mass[rm]:.3f}\n")
        
        # Atoms section
        file.write("\nAtoms\n\n")
        for index, row in atoms.iterrows():
            file.write(f"{int(row['Atom Number'])} {int(row['Chain'])} {int(row['Amino Acid'])} {int(row['Ignore'])} {row['x']:.3f} {row['y']:.3f} {row['z']:.3f}\n")

        # Bonds section
        file.write("\nBonds\n\n")
        for bond in bonds:
            bond_number, bond_type, atom1, atom2 = bond
            file.write(f"{bond_number} {bond_type} {atom1} {atom2}\n")

def main():
    """Main function to set up and run the simulation."""
    
    # Parse command-line arguments for temperature and model
    parser = argparse.ArgumentParser(description="Run protein simulation with specified parameters.")
    parser.add_argument("--temperature", type=float, default=300, help="Temperature for the simulation (in Kelvin).")
    parser.add_argument("--protein_num", type=float, default=0, help="Index of FASTA sequeunce in file proteins.csv")
    parser.add_argument("--model", type=str, default='calvados', help="Model which should be used")
    args = parser.parse_args()

    filename = 'proteins.csv'
    Nchains = 200
    sigma = 10
    lxy = 15.0  # nm
    lz = 150.0  # nm
    sampling = 5e5

    if args.model == 'calvados':
        run_calvados(args.temperature,args.protein_num,filename,Nchains,sigma,lxy,lz,sampling) 
    elif args.model == 'mpipi':
        run_mpipi(args.temperature,args.protein_num,filename,Nchains,sigma,lxy,lz,sampling)

    


if __name__ == "__main__":
    main()
