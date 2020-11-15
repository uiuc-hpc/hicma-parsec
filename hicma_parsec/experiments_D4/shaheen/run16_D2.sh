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
#SBATCH --nodes=16
#SBATCH --time=1:00:00

typ=O3
#typ=release

export MPICH_MAX_THREAD_SAFETY=multiple
#export LD_LIBRARY_PATH=/lustre/project/k1205/lei/parsec/build-$typ/install/dplasma/lib/:$LD_LIBRARY_PATH

numnodes=16; nummpi=16; ntasks_per_node=1; nthread=32; p=4; q=4
core=31

NO=1  #No. of runs, NO

maxrank=50
genmaxrank=150
compmaxrank=150

home=/lustre/project/k1205/lei/lorapo

sruncmd="srun --job-name=lorapo-$_m-$_mb-$SLURM_JOB_NUM_NODES --hint=nomultithread \
            --nodes=$numnodes \
            --ntasks=$nummpi \
            --ntasks-per-node=$ntasks_per_node  --cpus-per-task=${nthread} --hint=nomultithread "

N=1080000
NB=3600

# First run
for subtile in 300; do
    for D in 2; do
        for band in 1;  do
            for loop in `seq 1 $NO`; do
                cmd="$sruncmd \
                     numactl --interleave=all \
                     $home/testing_dpotrf -P $p -Q $q -N $N -t $NB -e 1e-8 --maxrank $maxrank -j $N -v -c $core --genmaxrank $genmaxrank --compmaxrank $compmaxrank -D $D --band $band -F 0 -Y 1 -z $subtile"
                echo $cmd
                eval $cmd
            done
        done
    done
done
