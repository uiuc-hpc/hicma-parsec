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
NB_max=10800     # max NB in loop

filename=shaheen_D4_performance_nodes"$nodes"

home=/lustre/project/k1205/lei/lorapo

sruncmd="srun --job-name=lorapo-$_m-$_mb-$SLURM_JOB_NUM_NODES --hint=nomultithread \
            --nodes=$nodes \
            --ntasks=$nodes \
            --ntasks-per-node=$ntasks_per_node  --cpus-per-task=${nthread} --hint=nomultithread "

mb_opt=0

# First run
#for N in 1080000 `seq 2160000 2160000 6480000`; do
# Second run
for N in 8640000; do
    time_opt=1000000000
    # All 
    #for NB in 1200 1500 1800 2400 2700 3000 3375 3600 4000 4500 5000 5400 6000 7200 9000 10000 10800 12000 15000 18000 20000; do
    #for NB in 5000 5400 6000 7200 9000 10000 10800 12000 15000 18000 20000; do
    for NB in 4500; do 
        maxrank=$(($NB/2))
        genmaxrank=$(($NB/2))
        compmaxrank=$(($NB/2))

        if [ $NB -ge $mb_opt ]; then
            out=0.0

            echo ""
            echo "begin: mb_opt $mb_opt; time_opt $time_opt"

            for loop in `seq 0 $NO`; do
                cmd="$sruncmd \
                     numactl --interleave=all \
                     $home/testing_dpotrf -P $p -Q $q -N $N -t $NB -e 1e-8 --maxrank $maxrank -j $N -v -c $core --genmaxrank $genmaxrank --compmaxrank $compmaxrank -D 4 -Z 1 -F 0 -Y 1 -z 300 -E 1"
                echo $cmd
                eval $cmd 2>&1 | tee tmp_$nodes.txt

                # Check memory
		cat tmp_$nodes.txt | grep memory_for_matrix_allocation_per_node
		if [ $? -eq 0 ]; then
                    echo "N $N NB $NB: out_of_memory"
                    if [[ $time_opt != 1000000000 ]]; then
                        echo ""
                        echo "Time optimal: $N $time_opt $mb_opt $var out_of_memory"
                        echo "$N $time_opt $mb_opt $var out_of_memory" | tee -a result/$filename.csv
                    fi

		    exit 1 
		fi

                cat tmp_$nodes.txt >> result/shaheen_D4_performance_all.csv
                all_put=`cat tmp_$nodes.txt | grep R-LO`
                echo $all_put
                echo $all_put >> result/${filename}_all.csv
                out[$loop]=`echo $all_put | awk '{printf $42}'`
                echo "Time: $N $NB ${out[$loop]}"
             done

             # Sort
             for ii in `seq 0 $NO`; do
                 for jj in `seq 0 $NO`; do
                     if [ $(bc <<< "${out[$ii]} < ${out[$jj]}") -eq 1 ]; then
                         t=${out[$ii]}
                         out[$ii]=${out[$jj]}
                         out[$jj]=$t
                     fi
                 done
             done

            echo ""
            echo "sorted result"
            for ii in `seq 0 $NO`; do
               echo "${out[$ii]}"
            done
            echo ""

            sum=0.0
            for ii in `seq 0 $NO`; do
                sum=`echo "$sum + ${out[$ii]}" | bc`
            done
            avg=`echo "$sum/$(($NO+1))" | bc -l`

            echo "Time min: $N $NB min ${out[0]} max ${out[$NO]} avg $avg"
            if [[ ${out[0]} != "" ]]; then
                if [ $(bc <<< "${out[0]} < $time_opt") -eq 1 ]; then
                   time_opt=${out[0]}
                   time_max=${out[$NO]}
                   time_avg=$avg
                   mb_opt=$NB
                   var=`echo "($time_max-$time_opt)/$time_opt" | bc -l`
                else
                   echo ""
                   echo "Time optimal: $N $mb_opt opt $time_opt max $time_max avg $time_avg ; var $var"
                   echo "$N $time_opt $mb_opt $var" | tee -a result/$filename.csv
                   echo ""
                   break 
                fi

                if [ $NB -eq $NB_max ]; then
                   echo ""
                   echo "Time optimal: $N $time_opt $mb_opt $var not_reach_bottom"
                   echo "$N $time_opt $mb_opt $var not_reach_bottom" | tee -a result/$filename
                   echo ""
                fi

            fi

            if [[ $NB == $NB_max && $time_opt != 1000000000 ]]; then
                echo ""
                echo "Time optimal: $N $time_opt $mb_opt $var only_one_got"
                echo "$N $time_opt $mb_opt $var only_one_got" | tee -a result/$filename.csv
                echo ""
            fi

        fi
    done
done
