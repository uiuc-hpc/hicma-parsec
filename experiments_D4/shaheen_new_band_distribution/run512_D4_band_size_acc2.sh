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
#SBATCH --nodes=512
#SBATCH --time=24:00:00

typ=O3
#typ=release

export MPICH_MAX_THREAD_SAFETY=multiple
#export LD_LIBRARY_PATH=/lustre/project/k1205/lei/parsec/build-$typ/install/dplasma/lib/:$LD_LIBRARY_PATH

nodes=$SLURM_JOB_NUM_NODES
ntasks_per_node=1
nthread=32;
P=16
Q=32

home=/lustre/project/k1205/lei/lorapo

sruncmd="srun --job-name=lorapo-$SLURM_JOB_NUM_NODES --hint=nomultithread \
            --nodes=$nodes \
            --ntasks=$nodes \
            --ntasks-per-node=$ntasks_per_node  --cpus-per-task=${nthread} --hint=nomultithread"

maxrank=200
for N in 1080000 2160000; do

  if [ $N -eq 1080000 ]; then
    NB=1000
  else
    NB=1200
  fi

  file_name=shaheen_"$nodes"node_band_size_acc2.txt
  for band_size in `seq 1 17`; do
    cmd="$sruncmd \
         numactl --interleave=all \
         $home/testing_dpotrf -P $P -Q $Q -N $N -t $NB -e 1e-2 --maxrank $maxrank -j $N -v -c 31 --genmaxrank $maxrank --compmaxrank $maxrank -D 4 --band $band_size -F 0 -Y 1 -z 300 -E 0" 
    echo $cmd
    eval $cmd 2>&1 | tee -a result/$file_name
  done
done
