#!/bin/bash

# Slurm job options (name, compute nodes, job time)
#SBATCH --job-name=automaton
#SBATCH --time=0:30:0
#SBATCH --exclusive
#SBATCH --nodes=2



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




mpi_list=(1 4 8 12 16 20 24 32 36)

for process_num in ${mpi_list[@]}
do
  echo "Process = $process_num"
  srun --unbuffered --cpu-bind=core -n $process_num ./reactor config_simple.txt output_synchronize_no_openmp.txt 1234
done


