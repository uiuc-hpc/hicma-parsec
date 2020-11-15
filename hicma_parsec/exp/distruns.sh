#!/bin/bash
dry="-"
nprocs="16 32 64 128 256"
nprocs="64"
cases=
trace="-"
op="potrf"
que="workq"; timelimit="06:00:00"
if [ $# -eq 1 -o $# -eq 2 ]; then
    cases=$1
    if [ ! -f "$cases" ];then
        echo "Cases file does not exist: $cases" 
        exit
    else
        . $cases
    fi
    if [ $# -eq 2 ]; then
        dry=$2
    fi
else 
    echo "Usage: $0 casesFile [dry]" 1>&2
    exit
fi
#corespernode=40; numthreads=39 #skylake
corespernode=30; numthreads=31 #shaheen
#corespernode=30; numthreads="1 2 4 8 16 31" #shaheen

#que="debug"; timelimit="00:30:00"; echo "DEBUG DEBUG DEBUG DEBUG DEBUG "
for nodes in $nprocs; do
    echo "#Number of nodes: $nodes ============================="
    if [ "$dry" == "dry" ]; then
        cmdbatch=""
    else
        cmdbatch="sbatch --parsable --nodes=$nodes --partition=$que --time=$timelimit --account=k1205 --job-name=lo-$_appdata-$nodes \
           --output=/project/k1205/akbudak/lorapo/exp/out/%j \
           --error=/project/k1205/akbudak/lorapo/exp/err/%j \
           --cpus-per-task=32 \
           --mail-type=END,FAIL \
           --mail-user=kadir.akbudak@kaust.edu.sa \
            "
    fi
    for isub in 0;do
        sched="prio"

        prog="hic"
        caseids=${allcaseids[$nodes]}      
        ncases=`echo "$caseids" | wc -w`    
        startt=0;   endt=$((ncases-1)); 
        echo -n "#`date` on $nodes nodes. $note - $prog - $cases "
        echo -n \"$caseids\" 
        echo " $ncases ($startt-$endt-$step)"

        ct=$startt
        while [ $ct -le $endt ]; do
            et=$((ct+step-1))
            if [ $et -gt $endt ]; then
                et=$endt
            fi
            #ids="`seq $ct $et`"
            arrcaseids=($caseids)
            ids=
            for icase in `seq $ct $et`; do
                ids="$ids ${arrcaseids[icase]}"
            done
            ct=$((ct+step))
            echo "#case ids: $ids"
            for nt in $numthreads; do
                corespernode=$nt
                maxsub=$((nodes*(corespernode+10)))
                minsub=$((nodes*corespernode))
                $cmdbatch exp/distmem.sh $nodes $nodes $nt $trace $prog - "$ids" $dry $cases $timelimit $que $sched $maxsub $minsub $op
            done
        done
    done
done

