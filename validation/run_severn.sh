#!/bin/bash
#SBATCH -A prio_ey239
#SBATCH --job-name=severn_nest
#SBATCH --partition=priority_gpu
#SBATCH --gpus=h100:1
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=16G
#SBATCH --output=%x-%j.out
#SBATCH --mail-type=ALL
#SBATCH --exclude=a1127u05n01
nvidia-smi --query-gpu=index,uuid,memory.used,memory.total --format=csv
module purge
module load Julia

julia --project=test --threads=4 --gcthreads=4,1 validation/severn_estuary.jl ${1:-30}
