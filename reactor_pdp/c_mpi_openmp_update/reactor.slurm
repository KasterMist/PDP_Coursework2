#!/bin/bash

# Slurm job options (name, compute nodes, job time)
#SBATCH --job-name=reactor
#SBATCH --time=0:30:0
#SBATCH --exclusive
#SBATCH --nodes=12
#SBATCH --tasks-per-node=9
#SBATCH --cpus-per-task=4



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

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# mpi_list=()
# mpi_list=()
# mpi_list=(1 4 8 12 16 20 24 28 32 36 40 44 48 52 56 60 64 68 72 76 80 84 88 92 96 100 104 108)
# for process_num in ${mpi_list[@]}
# do
#   echo "Process = $process_num"
#   # srun --unbuffered --cpu-bind=core -n $process_num ./reactor config_simple.txt output_100syn-add.txt
#   srun -n $process_num ./reactor config_simple.txt output_1syn-add.txt 
# done

# echo "24"
srun --unbuffered --cpu-bind=core -n 108 ./reactor config_simple.txt output_weak_scaling.txt 1234



