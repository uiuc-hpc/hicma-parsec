include make.inc


#CC = mpicc
#LD = mpicc
PP = ${DPLASMA_INSTALL}/bin/parsec-ptgpp --noline
#PPFLAGS = -M index-array

ifdef DEBUG
CFLAGS += -g -O0
else
CFLAGS += -O3
endif

CFLAGS += -std=c11 -std=gnu99 -D_GNU_SOURCE
CFLAGS += -DGITHASH=$(shell git rev-parse HEAD)

# Add STARS-H lib
CFLAGS += $(shell PKG_CONFIG_PATH=$(PKG_CONFIG_PATH) pkg-config --cflags starsh)
LDFLAGS += $(shell PKG_CONFIG_PATH=$(PKG_CONFIG_PATH) pkg-config --libs starsh)

# Add shared libs
CFLAGS += -I$(DPLASMA_INSTALL)/include -I$(DPLASMA_INSTALL)/include/parsec -I${DPLASMA_SRC} -I${DPLASMA_SRC}/src/include -I${DPLASMA_SRC}/${BUILD_TYPE}/src/include -I${DPLASMA_SRC}/src -I${DPLASMA_SRC}/src/cores
LDFLAGS += -Wl,-rpath,$(DPLASMA_INSTALL)/lib64 $(DPLASMA_INSTALL)/lib64/libparsec.so
LDFLAGS += -Wl,-rpath,$(DPLASMA_INSTALL)/lib $(DPLASMA_INSTALL)/lib/libdplasma.so

# Add GPU
ifdef PARSEC_HAVE_CUDA 
LDFLAGS += -lcublas -lcudart -lcusolver
endif

# Add PaRSEC goodies
CFLAGS += $(shell PKG_CONFIG_PATH=$(DPLASMA_INSTALL)/lib/pkgconfig:$(PKG_CONFIG_PATH) pkg-config --cflags parsec)
LDFLAGS += $(shell PKG_CONFIG_PATH=$(DPLASMA_INSTALL)/lib/pkgconfig:$(PKG_CONFIG_PATH) pkg-config --libs parsec)
CFLAGS += -I${DPLASMA_SRC}/parsec -I${DPLASMA_SRC}/parsec/parsec
LDFLAGS += -L. -lhicmahcore

# Add DPLASMA goodies
CFLAGS += $(shell PKG_CONFIG_PATH=$(DPLASMA_INSTALL)/lib64/pkgconfig:$(DPLASMA_INSTALL)/lib/pkgconfig:$(PKG_CONFIG_PATH) pkg-config --cflags dplasma)
LDFLAGS += $(shell PKG_CONFIG_PATH=$(DPLASMA_INSTALL)/lib64/pkgconfig:$(DPLASMA_INSTALL)/lib/pkgconfig:$(PKG_CONFIG_PATH) pkg-config --libs dplasma) -L. -lhicmahcore
CFLAGS += -DPLASMA

JDFS += HiCMA_dpotrf_L_3flow.jdf HiCMA_dpotrf_L_2flow.jdf STARSH_gen.jdf STARSH_check.jdf matrix_gather.jdf band_free.jdf band_size_calculate.jdf band_gen.jdf rank_check.jdf reorder_gemm.jdf rank_gather.jdf Av_memory.jdf
GENERATED_SRC = $(subst .jdf,.c,${JDFS})
GENERATED_HDR = $(subst .jdf,.h,${JDFS})
GENERATED_OBJ = $(subst .jdf,.o,${JDFS})

.PRECIOUS: %.c %.o %.h #do not delete intermediate files

SRC += HiCMA_dpotrf_L_3flow_wrapper.c HiCMA_dpotrf_L_2flow_wrapper.c hicma_parsec.c testing_dpotrf.c
OBJS = $(subst .c,.o,${SRC})


all: libhicmahcore.a testing_dpotrf testing_tlr_gemm

libhicmahcore.a: hicma_hcore.c hicma_hcore.h hicma_init.c descprint.c
	$(CC) $(CFLAGS) -c hicma_hcore.c
	$(CC) $(CFLAGS) -c hicma_init.c
	$(CC) $(CFLAGS) -c descprint.c
	ar rc libhicmahcore.a hicma_hcore.o hicma_init.o descprint.o

%.c %.h: %.jdf Makefile
	$(PP) ${PPFLAGS} -E -i $< -o $(basename $<)

%.o: %.c ${GENERATED_HDR} Makefile
	$(CC) -c $< -o $@ $(CFLAGS)

testing_dpotrf: ${GENERATED_OBJ} ${OBJS}
	@echo
	@echo CFLAGS: $(CFLAGS)
	@echo
	@echo LDFLAGS: $(LDFLAGS)
	@echo
	$(LD) -o $@ $^ ${LDFLAGS}

testing_tlr_gemm.o: testing_tlr_gemm.c Makefile
	$(CC) -c $< -o $@ $(CFLAGS)

testing_tlr_gemm: testing_tlr_gemm.o
	$(LD) -o $@ $^ ${LDFLAGS}

clean:
	rm -f ${GENERATED_OBJ} ${GENERATED_HDR} ${GENERATED_SRC} ${OBJS} testing_dpotrf testing_tlr_gemm hicma_hcore.o libhicmahcore.a hicma_init.o descprint.o *.o
