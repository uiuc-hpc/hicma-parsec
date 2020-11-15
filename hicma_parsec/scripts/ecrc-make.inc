CC = mpicc
LD = mpicc
PARSEC_SRC = ../parsec
PARSEC_ROOT = /home/akbudak/lib/parsecnoplasma/
ifdef PLASMA
    PARSEC_ROOT = /home/akbudak/lib/parsec/
endif
PARSEC_BIN = ${PARSEC_ROOT}
DPLASMA_DIR = ${PARSEC_BIN}/dplasma
PLASMA_DIR=../plasma
