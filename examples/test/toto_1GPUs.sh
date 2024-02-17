#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --nodes=1 --ntasks-per-node=1 --cpus-per-task=32
###SBATCH --out=job_1GPUs.out
###SBATCH --gres=gpu:4
#SBATCH --gpus-per-node=1
#SBATCH --account=FUSIO_HLST_4
###SBATCH --account=FUAL7_TSVV3       # account name
#SBATCH --partition=boost_fua_prod # partition name
###SBATCH --partition=boost_fua_dbg # partition name
###SBATCH  -A FUSIO_HLST_4
###export CUDA_DEVICE_ORDER="PCI_BUS_ID"
export NV_ACC_NOTIFY=1
ulimit -s unlimited
#export NV_ACC_CUDA_STACKSIZE=64000000
#export CU_LIMIT_STACK_SIZE=256000000
mpirun --rank-by core ./mpirun_bind.sh  /leonardo_scratch/large/userexternal/mpeybern/shyLab/build/shy triple_point_multimat_LF.json
