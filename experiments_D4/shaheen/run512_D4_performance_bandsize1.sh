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

export MPICH_MAX_THREAD_SAFETY=multiple
ntasks_per_node=1
nthread=32

nodes=$SLURM_JOB_NUM_NODES

if [ $nodes -eq 16 ]; then
    p=4
    q=4
elif [ $nodes -eq 32 ]; then
    p=4
    q=8
elif [ $nodes -eq 64 ]; then
    p=8
    q=8
elif [ $nodes -eq 128 ]; then
    p=8
    q=16
elif [ $nodes -eq 256 ]; then
    p=16
    q=16
elif [ $nodes -eq 512 ]; then
    p=16
    q=32
elif [ $nodes -eq 1024 ]; then
    p=32
    q=32
else
    echo "Wrong node count"
fi

core=31
NO=0  #No. of runs, NO+1

filename=shaheen_D4_performance_nodes"$nodes"_bandsize1_2

home=/lustre/project/k1205/lei/lorapo

sruncmd="srun --job-name=lorapo-$_m-$_mb-$SLURM_JOB_NUM_NODES --hint=nomultithread \
            --nodes=$nodes \
            --ntasks=$nodes \
            --ntasks-per-node=$ntasks_per_node  --cpus-per-task=${nthread} --hint=nomultithread "

# First run
#for N in 1080000 2160000 3240000; do
for N in 4320000; do 
    if [ $N -eq 1080000 ]; then
        NB=2400
    elif [ $N -eq 2160000 ]; then
        NB=2700
    elif [ $N -eq 3240000 ]; then
        NB=3000
    else
        NB=3375
    fi

    maxrank=$(($NB/2))
    genmaxrank=$(($NB/2))
    compmaxrank=$(($NB/2))

    for loop in `seq 0 $NO`; do
        cmd="$sruncmd \
             numactl --interleave=all \
             $home/testing_dpotrf -P $p -Q $q -N $N -t $NB -e 1e-8 --maxrank $maxrank -j $N -v -c $core --genmaxrank $genmaxrank --compmaxrank $compmaxrank -D 4 -Z 1 -F 0 -Y 1 -z 300 -E 0"
        echo $cmd
        eval $cmd 2>&1 | tee -a result/$filename.txt
    done
done
