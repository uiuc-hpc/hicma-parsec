#_appdata="--ss"; timelimit="10:00:00"
_appdata="3" #sqexp 3d 
timelimit="16:00:00"
#timelimit="08:00:00"
timelimit="04:00:00"
#timelimit="00:30:00" #TODO
tlmxrk=200
_compmaxrank=200 
_genmaxrank=200
_appdata=a; 
_wavek=50; _genmaxrank=50; tlmxrk=50; _compmaxrank=50 
_wavek=100; _genmaxrank=100; tlmxrk=100; _compmaxrank=100 

mrwaves[50]=50
mrwaves[100]=100
mrapps[2]=50
mrapps[3]=200
mrapps[4]=1680
#heatmap
#_subtile=300
allcaseids[16]="1 2 3 4 5"
nprocs="16"
_ci=1 ;apps[$_ci]="1"; waveks[$_ci]="50";tlmxrk=${mrwaves[50]}; nrows[$_ci]=1080000;  nb[$_ci]=2700;    acc[$_ci]=8;    maxrank[$_ci]=$tlmxrk;  genmaxrank[$_ci]=$tlmxrk;  compmaxrank[$_ci]=$tlmxrk; bands[$_ci]="0"; send_full_tiles[$_ci]="0"; lookaheads[$_ci]="0"; subtiles[$_ci]="${nb[_ci]} $_subtile"
_ci=2 ;apps[$_ci]="1"; waveks[$_ci]="100";tlmxrk=${mrwaves[100]}; nrows[$_ci]=1080000;  nb[$_ci]=2700;    acc[$_ci]=8;    maxrank[$_ci]=$tlmxrk;  genmaxrank[$_ci]=$tlmxrk;  compmaxrank[$_ci]=$tlmxrk; bands[$_ci]="0"; send_full_tiles[$_ci]="0"; lookaheads[$_ci]="0"; subtiles[$_ci]="${nb[_ci]} $_subtile"
_ci=3 ;apps[$_ci]="2"; waveks[$_ci]="0";tlmxrk=${mrapps[2]};  nrows[$_ci]=1080000;  nb[$_ci]=2700;    acc[$_ci]=8;    maxrank[$_ci]=$tlmxrk;  genmaxrank[$_ci]=$tlmxrk;  compmaxrank[$_ci]=$tlmxrk; bands[$_ci]="0"; send_full_tiles[$_ci]="0"; lookaheads[$_ci]="0"; subtiles[$_ci]="${nb[_ci]} $_subtile"
_ci=4 ;apps[$_ci]="3"; waveks[$_ci]="0";tlmxrk=${mrapps[3]};  nrows[$_ci]=1080000;  nb[$_ci]=2700;    acc[$_ci]=8;    maxrank[$_ci]=$tlmxrk;  genmaxrank[$_ci]=$tlmxrk;  compmaxrank[$_ci]=$tlmxrk; bands[$_ci]="0"; send_full_tiles[$_ci]="0"; lookaheads[$_ci]="0"; subtiles[$_ci]="${nb[_ci]} $_subtile"
_ci=5 ;apps[$_ci]="4"; waveks[$_ci]="0";tlmxrk=${mrapps[4]};  nrows[$_ci]=1080000;  nb[$_ci]=2700;    acc[$_ci]=8;    maxrank[$_ci]=$tlmxrk;  genmaxrank[$_ci]=$tlmxrk;  compmaxrank[$_ci]=$tlmxrk; bands[$_ci]="0"; send_full_tiles[$_ci]="0"; lookaheads[$_ci]="0"; subtiles[$_ci]="${nb[_ci]} $_subtile"

#breakdown
_subtile=300
allcaseids[64]="`seq 1 15`"
allcaseids[128]="`seq 1 15`"
allcaseids[128]="10 15"
nprocs="64"
nprocs="128"
nbms[1080000]=2700
nbms[2295000]=3375
nbms[3510000]=4500
nbms[2160000]=2700
nbms[3240000]=5400
_ci=0
for _m in 1080000 2160000 3240000 ;do #  2295000 3510000;do 
    _nb=${nbms[_m]}
    _ci=$((_ci+1)) ;apps[$_ci]="1"; waveks[$_ci]="50";tlmxrk=${mrwaves[50]}; nrows[$_ci]=$_m;  nb[$_ci]=$_nb;    acc[$_ci]=8;    maxrank[$_ci]=$tlmxrk;  genmaxrank[$_ci]=$tlmxrk;  compmaxrank[$_ci]=$tlmxrk; bands[$_ci]="1"; send_full_tiles[$_ci]="0"; lookaheads[$_ci]="1"; subtiles[$_ci]="$_subtile"
    _ci=$((_ci+1)) ;apps[$_ci]="1"; waveks[$_ci]="100";tlmxrk=${mrwaves[100]}; nrows[$_ci]=$_m;  nb[$_ci]=$_nb;    acc[$_ci]=8;    maxrank[$_ci]=$tlmxrk;  genmaxrank[$_ci]=$tlmxrk;  compmaxrank[$_ci]=$tlmxrk; bands[$_ci]="1"; send_full_tiles[$_ci]="0"; lookaheads[$_ci]="1"; subtiles[$_ci]="$_subtile"
    _ci=$((_ci+1)) ;apps[$_ci]="2"; waveks[$_ci]="0";tlmxrk=${mrapps[2]};  nrows[$_ci]=$_m;  nb[$_ci]=$_nb;    acc[$_ci]=8;    maxrank[$_ci]=$tlmxrk;  genmaxrank[$_ci]=$tlmxrk;  compmaxrank[$_ci]=$tlmxrk; bands[$_ci]="1"; send_full_tiles[$_ci]="0"; lookaheads[$_ci]="1"; subtiles[$_ci]="$_subtile"
    _ci=$((_ci+1)) ;apps[$_ci]="3"; waveks[$_ci]="0";tlmxrk=${mrapps[3]};  nrows[$_ci]=$_m;  nb[$_ci]=$_nb;    acc[$_ci]=8;    maxrank[$_ci]=$tlmxrk;  genmaxrank[$_ci]=$tlmxrk;  compmaxrank[$_ci]=$tlmxrk; bands[$_ci]="1"; send_full_tiles[$_ci]="0"; lookaheads[$_ci]="1"; subtiles[$_ci]="$_subtile"
    _ci=$((_ci+1)) ;apps[$_ci]="4"; waveks[$_ci]="0";tlmxrk=${mrapps[4]};  nrows[$_ci]=$_m;  nb[$_ci]=$_nb;    acc[$_ci]=8;    maxrank[$_ci]=$tlmxrk;  genmaxrank[$_ci]=$tlmxrk;  compmaxrank[$_ci]=$tlmxrk; bands[$_ci]="1"; send_full_tiles[$_ci]="0"; lookaheads[$_ci]="1"; subtiles[$_ci]="$_subtile"
done
step=1
note="Hicma $_appdata - $sizes - wavek:$_wavek - timelimit:$timelimit - compmaxrank:$_compmaxrank - tile max rank:$tlmxrk - genmaxrank:$_genmaxrank "

