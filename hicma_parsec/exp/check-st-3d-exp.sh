#!/bin/bash
#SBATCH --job-name=lorapo
#SBATCH --account=k1205
#SBATCH --output=/project/k1205/akbudak/lorapo/exp/out/%j
#SBATCH --error=/project/k1205/akbudak/lorapo/exp/err/%j
#SBATCH --cpus-per-task=32
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kadir.akbudak@kaust.edu.sa
#SBATCH --threads-per-core=1
#SBATCH --hint=nomultithread
#SBATCH --partition=debug
#SBATCH --ntasks-per-node=1
#SBATCH --time=00:10:59
#SBATCH --nodes=4

. /lustre/project/k1205/akbudak/lorapo/scripts/modules-shaheen-intel.sh
#srun /lustre/project/k1205/akbudak/lorapo/testing_dpotrf -P 2 -Q 2 5000 -t 500 -e 1e-8 -w 1 -u 400 -j 5000
#exit

#numnodes=1; nummpi=1; ntasks_per_node=1; nthread=1; __p=1; __q=1
numnodes=4; nummpi=4; ntasks_per_node=1; nthread=32; __p=2; __q=2
#numnodes=16; nummpi=16; ntasks_per_node=1; nthread=32; __p=4; __q=4
sruncmd="srun --job-name=lorapo-$_m-$_mb-$SLURM_JOB_NUM_NODES --hint=nomultithread \
            --nodes=$numnodes \
            --ntasks=$nummpi \
            --ntasks-per-node=$ntasks_per_node  --cpus-per-task=${nthread} --hint=nomultithread "
if [ \
    "$HOSTNAME" == "almaha" \
    -o  "$HOSTNAME" == "kw60319" \
    ]; then
    sruncmd="mpirun -n $numnodes " 
    #sruncmd="gdb -ex run --args"
fi
exps="5"
i=0; problem[$i]=1; m[$i]=5000;  mb[$i]=500;  w[$i]=50; stmxrk[$i]=$((mb[$i]/2)); shmxrk[$i]=$((mb[$i]/2)); cpmxrk[$i]=$((mb[$i]/2))
i=1; problem[$i]=2; m[$i]=5000;  mb[$i]=500;  w[$i]=-1; stmxrk[$i]=$((mb[$i]/2)); shmxrk[$i]=$((mb[$i]/2)); cpmxrk[$i]=$((mb[$i]/2))
i=2; problem[$i]=3; m[$i]=10000; mb[$i]=2000; w[$i]=-1; stmxrk[$i]=$((mb[$i]/2)); shmxrk[$i]=$((mb[$i]/2)); cpmxrk[$i]=$((mb[$i]))
i=3; problem[$i]=4; m[$i]=10000; mb[$i]=2000; w[$i]=-1; stmxrk[$i]=$((mb[$i]/2)); shmxrk[$i]=$((mb[$i]/2)); cpmxrk[$i]=$((mb[$i]))
i=4; problem[$i]=5; m[$i]=10000; mb[$i]=2000; w[$i]=50; stmxrk[$i]=$((mb[$i]/2)); shmxrk[$i]=$((mb[$i]/2)); cpmxrk[$i]=$((mb[$i]/2))
i=5; problem[$i]=4; m[$i]=2500;  mb[$i]=500;  w[$i]=-1; stmxrk[$i]=$((mb[$i]/2));   shmxrk[$i]=$((mb[$i]));   cpmxrk[$i]=$((mb[$i]))


for iexp in $exps;do
    __problem=${problem[iexp]}
    __m=${m[iexp]}
    __mb=${mb[iexp]}
    __w=${w[iexp]}
    __stmxrk=${stmxrk[iexp]}
    __shmxrk=${shmxrk[iexp]}
    __cpmxrk=${cpmxrk[iexp]}
    cmd="$sruncmd \
        numactl --interleave=all \
        $PWD/testing_dpotrf -P $__p -Q $__q $__m \
        --cores $nthread \
        --MB $__mb \
        --fixedacc 1e-5 \
        --wavek $__w  \
        --adddiag $__m \
        --kind_of_problem $__problem \
        -F 0 \
        --band 2 \
        --maxrank $__stmxrk \
        --genmaxrank $__shmxrk \
        --compmaxrank $__cpmxrk \
        "
    echo "-------------------RUNNING WITHOUT ANY CHECK-------------------"
    echo $cmd
    #eval $cmd
    echo
    echo "-------------------RUNNING WITH CHECK-------------------"
    echo "$cmd --check"
    eval "$cmd --check"
done
#--check
