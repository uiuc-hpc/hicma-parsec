#!/bin/bash

numnodes=1
numtasks=1
numthreads=31
tracearg=""
prog=""
minmaxsub=""
dry=""
if [ $# -eq 15 ]; then
    numnodes=$1
    nummpi=$2
    numthreads=$3
    trace=$4
    prog=$5
    minmaxsub=$6
    exps="$7"
    dry=$8
    sizefile=$9
    time_limit=${10}
    queue=${11}
    scheds=${12}
    maxsub=${13}
    minsub=${14}
    op=${15}
else
    echo "usage: $0 numnodes nummpi numthreads [trace,-] [hic,cham] [enable minmaxsub depending on #tasks,-] exps [dry,-] sizefile time_limit queue scheds [maxsub,-] [minsub,] [potrf posv]"
    echo
    echo "Your input:"
    echo $*
    exit
fi
#Calculate number of nodes and number of mpi tasks per node
sr=$numnodes #$(echo "scale=0;sqrt ( $numnodes ) / 1" | bc -l) ;
p=$(echo " ( l($sr) / l(2) )/2" | bc -l) ;
p=$(echo " scale=0;( $p / 1 ) " | bc -l) ;
sqrt_numnodes=$(echo " scale=0;( 2 ^ $p )/1 " | bc -l) ;
sqrt_numnodesQ=$((numnodes/sqrt_numnodes))
#echo $sqrt_numnodes $sqrt_numnodesQ; exit
ntasks_per_node=$((nummpi/numnodes))

if [ "$op" == "potrf" ]; then
    BINDIR=/project/k1205/akbudak/lorapo
    #BINDIR=/lustre/project/k1205/peiy/lorapo/
    BINNAME=testing_d${op}
else
    echo "Wrong op: $op"
    exit
fi
BIN=$BINDIR/$BINNAME

if [ "$dry" != "dry" ]; then
    module load cdt
    . /lustre/project/k1205/akbudak/lorapo/scripts/modules-shaheen-intel.sh
fi
#echo "======================================================================="
#echo "||   DO NOT FORGET TO USE  cc -craympich-mt  WHILE COMPILING PARSEC  ||"
#echo "======================================================================="
export MPICH_MAX_THREAD_SAFETY=multiple
#export LD_LIBRARY_PATH=/lustre/project/k1205/akbudak/codes/parsec/build-intel/install/lib/:/lustre/project/k1205/akbudak/codes/parsec/build-intel/install/dplasma/lib/:$LD_LIBRARY_PATH
export PKG_CONFIG_PATH=/lustre/project/k1205/peiy/parsec/build/install/dplasma/lib/pkgconfig/:/lustre/project/k1205/peiy/parsec/build/install/lib/pkgconfig:/lustre/project/k1205/akbudak/codes/plasma-installer_2.8.0/install-intel/lib/pkgconfig:/lustre/project/k1205/akbudak/hicma-dev/stars-h/build/install/lib/pkgconfig:$PKG_CONFIG_PATH
export LD_LIBRARY_PATH=/lustre/project/k1205/peiy/parsec/build/install/lib/:/lustre/project/k1205/peiy/parsec/build/install/dplasma/lib:$LD_LIBRARY_PATH

#Print linked libraries to stderr
if [ "$dry" != "dry" ]; then
    ldd $BIN 1>&2
    echo $LD_LIBRARY_PATH  1>&2
fi

#Parameters for HiCMA
_factor=4       #ratio of mb/maxrank
_acc=8          #fixed accuracy


#Sizes of matrices for experimental cases
if [ ! -f "$sizefile" ]; then
    echo "Size file does not exist: $sizefile"
    exit
fi

#DEFAULT PARAMETERS
_appdata="--ss" #_maxrank=100 # for 54K
_appdata="--edsin";_wavek=40
_decay=0
_compmaxrank=1
#CUSTOM PARAMETERS
. $sizefile


#cham gives this error for prio:
#[starpu][_starpu_priority_push_task][assert failure] task priority 1464 is not between minimum -5 and maximum 5

sruncmd="srun \
--hint=nomultithread \
--nodes=$numnodes \
--ntasks=$nummpi \
--ntasks-per-node=$ntasks_per_node \
--threads-per-core=1 \
numactl --interleave=all \
" 
echo $sruncmd
echo
#Loop over scheduling algorithms of StarPU
for sched in $scheds;do
    export STARPU_SCHED=$sched
    #Loop over experimental cases
    for iexp in $exps;do 
        echo Experiment case:$iexp nrows:${nrows[iexp]} mb:${nb[iexp]}
        _m=${nrows[iexp]} 
        _b=${nb[iexp]}
        _apps=${apps[iexp]}
        _waveks=${waveks[iexp]}
        _maxrank=${maxrank[iexp]}
        _genmaxrank=${genmaxrank[iexp]}
        _compmaxrank=${compmaxrank[iexp]}
        _nrhs=${nrhs[iexp]}
        _bands=${bands[iexp]}
        _send_full_tiles=${send_full_tiles[iexp]}
        _lookaheads=${lookaheads[iexp]}
        _subtiles=${subtiles[iexp]}
        if [ ! -z "${acc[iexp]}" ]; then 
            _acc=${acc[iexp]}
        fi
        if [ ! -z "${decay[iexp]}" ]; then 
            _decay=${decay[iexp]}
        fi
        if [ ! -z "$_compmaxrank}" ]; then 
            #echo $_compmaxrank
            :
        elif [ ! -z "${compmaxrank[iexp]}" ]; then 
            _compmaxrank=${compmaxrank[iexp]}
        else
            ## calculate maxrank used for buffers during computation
            scaledb=$((_b/10))
            scaledmaxrank=$((_maxrank*4))
            val=$scaledmaxrank
            if [ $scaledb -lt $scaledmaxrank ]; then
                val=$scaledb
            fi
            if [ $val -le $_wavek ]; then
                val=$((_wavek+50))
            fi
            _compmaxrank=$val
            _compmaxrank=$((_b/2))
        fi
        if [ -z "$_m" ]; then 
            continue
        fi
        if [ -z "$_b" ]; then 
            continue
        fi
        _mb=$_b;
        _n=$((_m/_mb*_maxrank))
        _nb=$_maxrank

        if [ $sqrt_numnodes -eq $sqrt_numnodesQ ]; then
            pdims=$sqrt_numnodes 
        else
            pdims="$sqrt_numnodes $sqrt_numnodesQ" 
        fi
        for _appdata in $_apps; do
            for _wavek in $_waveks; do
                for pdim in $pdims; do
                    for band in $_bands; do
                        for send_full_tile in $_send_full_tiles; do
                            for lookahead in $_lookaheads; do
                                for subtile in $_subtiles; do
                                    if [ "$prog" == "hic" ]; then
                                        rankfile=/project/k1205/akbudak/hicma/exp/ranks/$prog-$sched-$_m-$_mb-$_nb-$numnodes-$nummpi-$numthreads-$SLURM_JOBID
                                        if [ "$op" == "potrf" ]; then
                                            cmd="$BIN \
                                                --m=$_m \
                                                --n_range=$_n:$_n \
                                                --k=$_m \
                                                --mb=$_mb \
                                                --nb=$_maxrank \
                                                --nowarmup \
                                                --threads=$numthreads \
                                                --p=$pdim \
                                                $tracearg \
                                                --rk=0 \
                                                --acc=$_acc \
                                                $_appdata \
                                                --starshwavek=$_wavek \
                                                --starshdecay=$_decay \
                                                --starshmaxrank=$_compmaxrank \
                                                --rankfile=$rankfile \
                                                "
                                            cmd="$BIN -P $pdim  $_m -t $_mb -e 1e-$_acc -w $_wavek -u $_maxrank -G $_genmaxrank -U $_compmaxrank -j $_m -D $_appdata -F $send_full_tile -Z $band -Y $lookahead -z $subtile -c $numthreads"
                                        fi
                                    fi

                                    msg="M:$_m N:$_n MB:$_mb NB:$_nb MAXRANK:$_maxrank DATE:`date` CMD:$cmd CASE:$sizefile" 
                                    echo "!BEGIN:" $msg 
                                    if [ "$dry" == "dry" ]; then
                                        #echo $cmd
                                        echo $sruncmd $cmd
                                    else
                                        echo "!BEGIN:" $msg 1>&2 
                                        tstart=$SECONDS
                                        $sruncmd $cmd
                                        tend=$SECONDS
                                        time_sec=$((tend-tstart))
                                        time_min=$((time_sec/60))
                                        time_hour=$((time_min/60))
                                        echo
                                        echo "!END:" $msg SECOND:$time_sec MINUTE:$time_min HOUR:$time_hour
                                        echo "!END:" $msg 1>&2 
                                    fi
                                    date
                                    echo
                                done
                            done
                        done
                    done
                done
            done
        done
    done
done

exit 0
    --printindex \
    --check \
    --printmat \
    --printindex \
    --trace \
--tag-output --timestamp-output 

$SLURM_JOBID    Job ID  5741192 $PBS_JOBID
$SLURM_JOB_NAME Job Name    myjob   $PBS_JOBNAME
$SLURM_SUBMIT_DIR   Submit Directory    /lustre/payerle/work    $PBS_O_WORKDIR
$SLURM_JOB_NODELIST Nodes assigned to job   compute-b24-[1-3,5-9],compute-b25-[1,4,8]   cat $PBS_NODEFILE
$SLURM_SUBMIT_HOST  Host submitted from login-1.deepthought2.umd.edu    $PBS_O_HOST
$SLURM_JOB_NUM_NODES    Number of nodes allocated to job    2   $PBS_NUM_NODES
$SLURM_CPUS_ON_NODE Number of cores/node    8,3 $PBS_NUM_PPN
$SLURM_NTASKS   Total number of cores for job???    11  $PBS_NP
$SLURM_NODEID   Index to node running on
relative to nodes assigned to job   0   $PBS_O_NODENUM
$PBS_O_VNODENUM Index to core running on
within node 4   $SLURM_LOCALID
$SLURM_PROCID   Index to task relative to job   0   $PBS_O_TASKNUM - 1

