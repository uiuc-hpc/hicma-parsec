#!/bin/bash
#SBATCH --job-name=lorapo
#SBATCH --account=k1205
#SBATCH --output=/project/k1205/lei/lorapo/experiments/%j
#SBATCH --error=/project/k1205/lei/lorapo/experiments/%j
#SBATCH --cpus-per-task=32
#SBATCH --threads-per-core=1
#SBATCH --hint=nomultithread
#SBATCH --partition=workq
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1024
#SBATCH --time=24:00:00

export MPICH_MAX_THREAD_SAFETY=multiple
export LD_LIBRARY_PATH=/lustre/project/k1205/lei/parsec/build-release/install/lib/:/lustre/project/k1205/lei/parsec/build-release/install/dplasma/lib/:$LD_LIBRARY_PATH

numnodes=1024; nummpi=1024; ntasks_per_node=1; nthread=32; p=32; q=32
core=31
Y=1   # lookahead
F=0   # not send full
band=1
subtile=300

D=3
maxrank=200
genmaxrank=300
compmaxrank=300

home=/lustre/project/k1205/lei/lorapo

sruncmd="srun --job-name=lorapo-$_m-$_mb-$SLURM_JOB_NUM_NODES --hint=nomultithread \
            --nodes=$numnodes \
            --ntasks=$nummpi \
            --ntasks-per-node=$ntasks_per_node  --cpus-per-task=${nthread} --hint=nomultithread "

# First run
for N in 24883200; do
    for NB in 10800; do
        memory_node=$(($numnodes*1024*1024*1024))
        NT=`echo "$N/$NB" | bc -l`
        Ag=`echo "$NT*$NB*$NB" | bc -l`
        uv=`echo "$NT*($NT-1)*$NB*$maxrank" | bc -l`
        memory=`echo "($Ag+$uv)/$memory_node*8" | bc -l`
        echo "N $N NB $NB needs_memory $memory G"

                cmd="$sruncmd \
                     numactl --interleave=all \
                     $home/testing_dpotrf -P $p -Q $q -N $N -t $NB -e 1e-8 --maxrank $maxrank -j $N --band $band -v -c $core --genmaxrank $genmaxrank --compmaxrank $compmaxrank -F $F -D $D -Y $Y -z $subtile"
                echo $cmd
                eval $cmd
    done
done
