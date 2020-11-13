#!/bin/bash
#SBATCH --job-name=lorapo
#SBATCH --account=k1339
#SBATCH --output=out/%j
#SBATCH --error=out/%j
#SBATCH --cpus-per-task=32
#SBATCH --threads-per-core=1
#SBATCH --hint=nomultithread
#SBATCH --partition=workq
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=128
#SBATCH --time=24:00:00

typ=O3
#typ=release

export MPICH_MAX_THREAD_SAFETY=multiple
#export LD_LIBRARY_PATH=/lustre/project/k1205/lei/parsec/build-$typ/install/dplasma/lib/:$LD_LIBRARY_PATH

nodes=$SLURM_JOB_NUM_NODES
ntasks_per_node=1
nthread=32;

home=/lustre/project/k1205/lei/lorapo

sruncmd="srun --job-name=lorapo-$SLURM_JOB_NUM_NODES --hint=nomultithread \
            --nodes=$nodes \
            --ntasks=$nodes \
            --ntasks-per-node=$ntasks_per_node  --cpus-per-task=${nthread} --hint=nomultithread"

N=2160000
P=8
Q=16
file_name=shaheen_"$nodes"node_"$N"_tile_size_new.txt

#for NB in 1800 2400 2700 3000 3375 3600 4000 4500 5000 5400 6000; do 
for NB in 1200 1000; do 
    maxrank=$(($NB/2))
    cmd="$sruncmd \
         numactl --interleave=all \
         $home/testing_dpotrf -P $P -Q $Q -N $N -t $NB -e 1e-8 --maxrank $maxrank -j $N -v -c 31 --genmaxrank $maxrank --compmaxrank $maxrank -D 4 --band 7 -F 0 -Y 1 -z 300 -E 1" 
    echo $cmd
    eval $cmd 2>&1 | tee -a result/$file_name
done
