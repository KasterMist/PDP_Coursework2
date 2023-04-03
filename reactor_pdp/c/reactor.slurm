#!/bin/bash

# Slurm job options (name, compute nodes, job time)
#SBATCH --job-name=automaton
#SBATCH --time=0:30:0
#SBATCH --exclusive
#SBATCH --nodes=1



#SBATCH --output=%x-%j.out

# Replace [budget] below with your unique budget code 
# For EPCC Students: [budget]=dc134-[student-id]
# For External Students: [budget]=dc135-[student-id]
#SBATCH --account=m22oc-s2312252

# We use the "standard" partition as we are running on CPU nodes
#SBATCH --partition=standard

# We use the "standard" QoS as our runtime is less than 4 days
#SBATCH --qos=standard

# Load any required modules
module load mpt
module load intel-compilers-19

# Change to the submission directory
cd $SLURM_SUBMIT_DIR





srun --unbuffered --cpu-bind=core -n 2 ./reactor config_simple.txt output.txt 1234


