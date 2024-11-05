# Protein Phase Diagram Simulation Pipeline

This project provides a streamlined, automated pipeline for generating phase diagrams of proteins using state-of-the-art coarse-grained molecular dynamics models. It supports simulations with the **Calvados** force field in OpenMM and the **mpipi** model in LAMMPS. The pipeline automates data preparation, simulation execution, checkpointing, and restarting, with flexible configuration for parameters like temperature, sampling, and protein sequences.

## Key Features

- **Dual Force Field Support**: Choose between:
  - **Calvados** in OpenMM: A force field tailored for protein coarse-graining.
  - **mpipi** in LAMMPS: A simplified model for protein simulations.
- **Automated Data Preparation**: Generates initial configurations with hexagonal distribution of protein chains, customizable chain count, box dimensions, and separation.
- **Checkpoint and Restart**: Save and resume simulations from checkpoints (`restart.equil` for mpipi; OpenMM checkpoints for Calvados).
- **Parameterizable**: Set simulation parameters, like temperature and sampling intervals, directly from the command line.
- **Flexible Simulation Management**: Includes PBS job submission scripts for cluster environments, enabling high-throughput simulation across parameter spaces.

## Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/yourusername/protein-phase-diagram-simulation.git
   cd protein-phase-diagram-simulation
