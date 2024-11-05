#!/bin/bash

# Define temperature range and step
start_temp=230
end_temp=330
num_temps=10
protein_num=2

# Define forcefields
forcefields=("calvados" "mpipi")

# Calculate temperature step size
temp_step=$(( (end_temp - start_temp) / (num_temps - 1) ))

# Create directories and run simulations
for forcefield in "${forcefields[@]}"; do
    # Create the main directory for each forcefield
    mkdir -p "$forcefield"

    # Loop over each temperature
    for i in $(seq 0 $((num_temps - 1))); do
        # Calculate the current temperature
        temp=$((start_temp + i * temp_step))

        # Create a directory for the specific temperature
        temp_dir="$forcefield/$temp"
        mkdir -p "$temp_dir"

        # Copy all necessary files to the temperature directory
        if [ "$forcefield" == "calvados" ]; then
            cp submit_calvados.pbs "$temp_dir/"
            cp calvados.py "$temp_dir/"
            cp mpipi_GG.py "$temp_dir/"
            cp proteins.csv "$temp_dir/"
            cp combined.py "$temp_dir/"

            # Navigate to the temperature directory
            cd "$temp_dir" || exit

            # Submit the job with temperature and forcefield model as arguments
            qsub -v TEMPERATURE="$temp",MODEL="$forcefield",PROTEIN_NUM="$protein_num" submit_calvados.pbs

            # Return to the root directory
            cd - > /dev/null

        elif [ "$forcefield" == "mpipi" ]; then
            cp submit_lammps.pbs "$temp_dir/"
            cp combined.py "$temp_dir/"
            cp mpipi_GG.py "$temp_dir/"
            cp calvados.py "$temp_dir/"
            cp proteins.csv "$temp_dir/"
            cp pairs.txt "$temp_dir/"

            # Navigate to the temperature directory
            cd "$temp_dir" || exit


            # Submit the job with temperature and forcefield model as arguments
            qsub -v TEMPERATURE="$temp",MODEL="$forcefield",PROTEIN_NUM="$protein_num" submit_lammps.pbs

            # Return to the root directory
            cd - > /dev/null
        fi
    done
done

echo "Simulation setup, configuration, and submission completed for all temperatures and forcefields."
