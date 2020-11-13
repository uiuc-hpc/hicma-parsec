	    CC = mpicc
	    LD = mpicc
            #NVCC = nvcc

	    BUILD_TYPE = build-O0

	    LDFLAGS = -lmkl_intel_lp64 -lmkl_sequential -lgsl

	    DPLASMA_SRC = /home/qcao3/lab/dplasma
            DPLASMA_BUILD = ${DPLASMA_SRC}/${BUILD_TYPE}
	    PARSEC_SRC = ${DPLASMA_SRC}/parsec
	    DPLASMA_INSTALL = ${DPLASMA_SRC}/${BUILD_TYPE}/install
