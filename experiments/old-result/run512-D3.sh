#!/bin/bash
#SBATCH --job-name=lorapo
#SBATCH --account=k1205
#SBATCH --output=/project/k1205/lei/lorapo/exp-utk/out/%j
#SBATCH --error=/project/k1205/lei/lorapo/exp-utk/out/%j
#SBATCH --cpus-per-task=32
#SBATCH --threads-per-core=1
#SBATCH --hint=nomultithread
#SBATCH --partition=workq
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=512
#SBATCH --time=24:00:00

export MPICH_MAX_THREAD_SAFETY=multiple
export LD_LIBRARY_PATH=/lustre/project/k1205/lei/parsec/build-release/install/lib/:/lustre/project/k1205/lei/parsec/build-release/install/dplasma/lib/:$LD_LIBRARY_PATH

numnodes=512; nummpi=512; ntasks_per_node=1; nthread=32; p=32; q=16
core=31
NO=0  #No. of runs, NO+1
Y=1   # lookahead
F=0   # not send full
band=1
subtile=300

D=3
maxrank=400
genmaxrank=500
compmaxrank=500

home=/lustre/project/k1205/lei/lorapo

sruncmd="srun --job-name=lorapo-$_m-$_mb-$SLURM_JOB_NUM_NODES --hint=nomultithread \
            --nodes=$numnodes \
            --ntasks=$nummpi \
            --ntasks-per-node=$ntasks_per_node  --cpus-per-task=${nthread} --hint=nomultithread "

NB_max=21600     # max NB in loop
mb_opt=0
memory_thd=250.0
# First run
for N in `seq 3240000 2160000 19440000`; do
    time_opt=1000000000
    # All 
    #for NB in 1200 1500 1800 2400 2700 3000 3375 3600 4000 4500 5000 5400 6000 7200 9000 10000 10800 12000 15000 18000 20000; do
    for NB in 2700 3000 3375 3600 4500 5000 5400 6000 7200 9000 10000 10800 12000 15000 18000 20000 21600; do
        if [ $NB -ge $mb_opt ]; then
            out=0.0
            memory_node=$(($numnodes*1024*1024*1024))
            NT=`echo "$N/$NB" | bc -l`
            Ag=`echo "$NT*$NB*$NB" | bc -l`
            uv=`echo "$NT*($NT-1)*$NB*$maxrank" | bc -l`
            memory=`echo "($Ag+$uv)/$memory_node*8" | bc -l`
            echo "N $N NB $NB needs_memory $memory G"

            if [ $(bc <<< "$memory_thd < $memory") -eq 1 ]; then
                echo "N $N NB $NB: out_of_memory"
                if [[ $time_opt != 1000000000 ]]; then
                    echo ""
                    echo "Time optimal: $N $time_opt $mb_opt $var out_of_memory"
                    echo "$N $time_opt $mb_opt $var out_of_memory" >> result/"node$numnodes"_"D$D"_"band$band"
                    #echo "$N $time_max $mb_opt $var out_of_memory" >> result/"node$numnodes"_"D$D"_"w$w"_max
                    #echo "$N $time_avg $mb_opt $var out_of_memory" >> result/"node$numnodes"_"D$D"_"w$w"_avg
                fi
                continue
            fi

            echo ""
            echo "begin: mb_opt $mb_opt; time_opt $time_opt"

            for loop in `seq 0 $NO`; do
                cmd="$sruncmd \
                     numactl --interleave=all \
                     $home/testing_dpotrf -P $p -Q $q -N $N -t $NB -e 1e-8 --maxrank $maxrank -j $N --band $band -v -c $core --genmaxrank $genmaxrank --compmaxrank $compmaxrank -F $F -D $D -Y $Y -z $subtile"
                echo $cmd
                all_put=`eval $cmd | grep R-LO`
                echo $all_put
                echo $all_put >> all_results_3D
                out[$loop]=`echo $all_put | awk '{printf $28}'`
                echo "$N ${out[$loop]} $NB $var" >> result/"node$numnodes"_"D$D"_"band$band"_all
                echo "Time: $N ${out[$loop]} $NB"
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
                   echo "$N $time_opt $mb_opt $var" >> result/"node$numnodes"_"D$D"_"band$band"
                   #echo "$N $time_max $mb_opt $var" >> result/"node$numnodes"_"D$D"_"w$w"_max
                   #echo "$N $time_avg $mb_opt $var" >> result/"node$numnodes"_"D$D"_"w$w"_avg
                   echo ""
                   break 
                fi

                if [ $NB -eq $NB_max ]; then
                   echo ""
                   echo "Time optimal: $N $time_opt $mb_opt $var not_reach_bottom"
                   echo "$N $time_opt $mb_opt $var not_reach_bottom" >> result/"node$numnodes"_"D$D"_"band$band"
                   #echo "$N $time_max $mb_opt $var not_reach_bottom" >> result/"node$numnodes"_"D$D"_"w$w"_max
                   #echo "$N $time_avg $mb_opt $var not_reach_bottom" >> result/"node$numnodes"_"D$D"_"w$w"_avg
                   echo ""
                fi

            fi

            if [[ $NB == $NB_max && $time_opt != 1000000000 ]]; then
                echo ""
                echo "Time optimal: $N $time_opt $mb_opt $var only_one_got"
                echo "$N $time_opt $mb_opt $var only_one_got" >> result/"node$numnodes"_"D$D"_"band$band"
                #echo "$N $time_max $mb_opt $var only_one_got" >> result/"node$numnodes"_"D$D"_"w$w"_max
                #echo "$N $time_avg $mb_opt $var only_one_got" >> result/"node$numnodes"_"D$D"_"w$w"_avg
                echo ""
            fi

        fi

    done
done
