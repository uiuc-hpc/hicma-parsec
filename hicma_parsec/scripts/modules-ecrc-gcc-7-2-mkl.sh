module load mkl/2018-update-2
module load ecrc-extras
module load cmake/3.11.1
module load gcc/7.2.0 openmpi/3.0.0-gcc-7.2.0
#. /opt/share/intel/2018/compilers_and_libraries/linux/mkl/bin/mklvars.sh intel64
module load hwloc/1.11.8-gcc-7.2.0
module load plasma/2.8.0-gcc-7.2.0-mkl


INSTALLDIR=/home/akbudak/lib/parsec_plasma_gcc_mkl
if [ ! -d $INSTALLDIR  ]; then
    echo -e "use\n\n\tscripts/parsec-install-ecrc-with-plasma-gcc-mkl.sh\n\nto install parsec on ECRC servers."
    return
fi
export PARSEC_ROOT=$INSTALLDIR
if [ -d $INSTALLDIR/bin ]; then
    export PATH=$INSTALLDIR/bin:$PATH
fi
if [ -d $INSTALLDIR/lib/pkgconfig ]; then
    export PKG_CONFIG_PATH=$INSTALLDIR/lib/pkgconfig:$PKG_CONFIG_PATH
fi
