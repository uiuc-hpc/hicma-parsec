/**
 * @copyright (c) 2020 King Abdullah University of Science and Technology (KAUST).
 *                     The Universiy of Tennessee and The Universiy of Tennessee Research Foundation.
 *                     All rights reserved.
 **/
#ifndef _HICMA_PARSEC_INTERNAL_H
#define _HICMA_PARSEC_INTERNAL_H 

/* PaRSEC hearders */
#include "parsec/data_dist/matrix/sym_two_dim_rectangle_cyclic.h"
#include "parsec/data_dist/matrix/sym_two_dim_rectangle_cyclic_band.h"
#include "parsec/data_dist/matrix/two_dim_rectangle_cyclic.h"
#include "parsec/data_dist/matrix/two_dim_rectangle_cyclic_band.h"
#include "parsec/execution_stream.h"
#include "parsec/runtime.h"
#include "parsec/profiling.h"
#include "parsec/parsec_internal.h"
#include "parsec/utils/debug.h"
#include "parsec/utils/zone_malloc.h"
#include "parsec/data_dist/matrix/matrix.h"
#include "parsec/private_mempool.h"
#include "parsec/data_internal.h"
#include "parsec/utils/mca_param.h"
#include "parsec/tests/interfaces/superscalar/common_timing.h"

/* DPLASMA headers */
#include "dplasma.h"
#include "include/dplasmajdf.h"
#include "include/dplasmaaux.h"
#include "tests/flops.h"

/* BLAS and LAPACKE headers */
#include <cblas.h>
#include <lapacke.h>
#include <core_blas.h>
#include <core_dblas.h>

/* Starsh headers */ 
#include <starsh-randtlr.h>
#include <starsh-electrodynamics.h>
#include <starsh-spatial.h>

/* Local headers */
#include "auxdescutil.h" // _printmat
#include "ka_counts.h"
#include "hicma_hcore.h"

/* System headers */
#include <time.h>
#include <sys/time.h>
#include <string.h>
#include <math.h>
#include <starsh.h>
#include <malloc.h>
#include <stdlib.h>
#include <stdio.h>

#ifdef PARSEC_HAVE_MPI
#include <mpi.h>
#endif

/* Recursive Kernel */
#define PARSEC_HAVE_RECURSIVE 1

/* Recursive headers */
#if defined(PARSEC_HAVE_RECURSIVE)
#include "parsec/data_dist/matrix/subtile.h"
#include "parsec/recursive.h"
#endif

#if defined(PARSEC_HAVE_CUDA)
#include "parsec/include/parsec/execution_stream.h"
#include "parsec/mca/device/cuda/device_cuda_internal.h"
#include "parsec/mca/device/cuda/device_cuda.h"
#include <cublas.h>
#define CUBLAS_V2_H_
#include <cusolverDn.h>
#include "parsec/utils/zone_malloc.h"
#endif

/* string print */ 
#define xstr(a) #a 

/* Fluctuation, used in auto-band tuning */
#define FLUCTUATION 0.66666667

/* Show rank statistics */
#define PRINT_RANK 0

/* Statistic rank info in dcRank
 * Need PRINT_RANK 1
 * [0] : the initial rank distribution
 * [1] : the minimum rank distribution during factorization
 * [2] : the maximum rank distribution during factorization
 * [3] : the final rank distribution
 */
#define RANK_MAP_BUFF 4

/* Statistic rank info in dcRank
 * Need PRINT_RANK 1
 * [0] : the initial rank distribution
 * [1] : the minimum rank distribution during factorization
 * [2] : the maximum rank distribution during factorization
 * [3] : the final rank distribution
 * [4] : [2] - [1] 
 * [5] : All, 0 - 4 
 */
#define RANK_MAP_TYPE 0

/* Print thread execution time */
#define PRINT_THREAD_EXE_TIME 0

/* Print starting time of task in the critical path */ 
#define PRINT_CRITICAL_PATH_TIME 0

/* Contiguous memory for band
 * Maybe benefit for GPU, as pin memory between GPU and CPU
 */ 
#define BAND_MEMORY_CONTIGUOUS 0

/* Memory threshold */
#define THRESHOLD_MEMORY_PER_NODE 70

/* Using MPI routines to gather rank
 * only to calculate the average, min, and max rank,
 * because the rank order for each tile is not correct.
 * Otherwise, use JDF to gather, but much slover
 * when matrix is large
 */
#define GATHER 1

/* Print more info for debugging */
#define DEBUG_INFO 0

/* Update PASTE_CODE_PROGRESS_KERNEL below if you change this list */
enum iparam_t {
  IPARAM_RANK,         /* Rank                              */
  IPARAM_NNODES,       /* Number of nodes                   */
  IPARAM_NCORES,       /* Number of cores                   */
  IPARAM_NGPUS,        /* Number of GPUs                    */
  IPARAM_P,            /* Rows in the process grid          */
  IPARAM_Q,            /* Columns in the process grid       */
  IPARAM_N,            /* Number of columns of the matrix   */
  IPARAM_NB,           /* Number of columns in a tile       */
  IPARAM_HMB,          /* Small MB for recursive hdags */
  IPARAM_HNB,          /* Small NB for recursive hdags */
  IPARAM_CHECK,        /* Checking activated or not         */
  IPARAM_CHECKINV,     /* Inverse Checking activated or not */
  IPARAM_ASYNC,        /* Bench the asynchronous version    */
  IPARAM_VERBOSE,      /* How much noise do we want?        */
  IPARAM_LOWLVL_TREE,  /* Tree used for reduction inside nodes  (specific to xgeqrf_param) */
  IPARAM_HIGHLVL_TREE, /* Tree used for reduction between nodes (specific to xgeqrf_param) */
  IPARAM_QR_TS_SZE,    /* Size of TS domain                     (specific to xgeqrf_param) */
  IPARAM_QR_HLVL_SZE,  /* Size of the high level tree           (specific to xgeqrf_param) */
  IPARAM_RANDOM_SEED,  /* Seed for the pseudo-random generators */
  IPARAM_MATRIX_INIT,  /* Matrix generator type */
  IPARAM_QR_DOMINO,    /* Enable/disable the domino between the upper and the lower tree (specific to xgeqrf_param) */
  IPARAM_QR_TSRR,      /* Enable/disable the round-robin on TS domain */
  IPARAM_BUT_LEVEL,    /* Butterfly level */
  IPARAM_SCHEDULER,    /* User-selected scheduler */
  /* HiCMA options */
  IPARAM_FIXED_RANK,
  IPARAM_MAX_RANK,
  IPARAM_GEN_MAX_RANK,
  IPARAM_COMP_MAX_RANK,
  IPARAM_BAND,         /* Use band distribution */
  IPARAM_KIND_OF_PROBLEM,
  IPARAM_SEND_FULL_TILE,
  IPARAM_LOOKAHEAD,
  IPARAM_AUTO_BAND,
  IPARAM_REORDER_GEMM,
  IPARAM_TWO_FLOW,
  IPARAM_SIZEOF
};

enum dparam_t {
  DPARAM_FIXED_ACC,
  DPARAM_WAVEK,
  DPARAM_ADD_DIAG,
  DPARAM_SIZEOF
};

#define PASTE_CODE_IPARAM_LOCALS(iparam)                                \
    int rank  = iparam[IPARAM_RANK];                                    \
    int nodes = iparam[IPARAM_NNODES];                                  \
    int cores = iparam[IPARAM_NCORES];                                  \
    int gpus  = iparam[IPARAM_NGPUS];                                   \
    int P     = iparam[IPARAM_P];                                       \
    int Q     = iparam[IPARAM_Q];                                       \
    int N     = iparam[IPARAM_N];                                       \
    int NB    = iparam[IPARAM_NB];                                      \
    int HNB   = iparam[IPARAM_HNB];                                     \
    int NT    = (N%NB==0) ? (N/NB) : (N/NB+1);                          \
    int check = iparam[IPARAM_CHECK];                                   \
    int loud  = iparam[IPARAM_VERBOSE];                                 \
    double wave_k = dparam[DPARAM_WAVEK];                              \
    int maxrank = iparam[IPARAM_MAX_RANK];                             \
    int genmaxrank = iparam[IPARAM_GEN_MAX_RANK];                      \
    int compmaxrank = iparam[IPARAM_COMP_MAX_RANK];                    \
    int fixedrk = iparam[IPARAM_FIXED_RANK];                           \
    int band_size = iparam[IPARAM_BAND];                               \
    int lookahead = iparam[IPARAM_LOOKAHEAD];                          \
    int kind_of_problem = iparam[IPARAM_KIND_OF_PROBLEM];              \
    int send_full_tile = iparam[IPARAM_SEND_FULL_TILE];                \
    int reorder_gemm = iparam[IPARAM_REORDER_GEMM];                    \
    int two_flow = iparam[IPARAM_TWO_FLOW];                            \
    int auto_band = iparam[IPARAM_AUTO_BAND];

#define PASTE_CODE_DPARAM_LOCALS(dparam) \
    double tol = dparam[DPARAM_FIXED_ACC];     \
    double add_diag = dparam[DPARAM_ADD_DIAG];

/* Define a double type which not pass through the precision generation process */
#define PASTE_CODE_FLOPS( FORMULA, PARAMS ) \
  double gflops = -1.0, flops = FORMULA PARAMS;

/**
 *  * No macro with the name max or min is acceptable as there is
 *   * no way to correctly define them without borderline effects.
 *    */
#undef max
#undef min
static inline int my_own_max(int a, int b) { return a > b ? a : b; }
static inline int my_own_min(int a, int b) { return a < b ? a : b; }

/* starsh parameter structure */
typedef struct starsh_params_s
{
    void *data;
    STARSH_kernel *kernel;
    STARSH_int *index;
} starsh_params_t;

/* 2flow version non-blocking 
 *
 * @param [in] parsec:             parsec context 
 * @param [inout] A:               the data, already distributed and allocated
 * @param [inout] Ar:              the rank info, already distributed and allocated
 * @param [in] RG:                 used for reordering GEMM 
 * @param [in] Rank:               used for gathering time info during factorization 
 * @param [in] acc:                accuracy threshold 
 * @param [in] fixed_rank:         fixed rank threshold used in recompression stage of HCORE_GEMM 
 * @param [in] maxrank:            max rank threshold for storage
 * @param [in] lookahead:          lookahead to guide workflow 
 * @param [in] band_size:          band size to control dense tiles; band_size=1: only diagonal in dense 
 * @param [in] hmb:                hierarchical mb; the sub tile size used in the recursive kernel 
 * @param [in] compmaxrank:        max rank threshold used in computation 
 * @param [in] send_full_tile:     whether send full tile during factorization; default is false to give better performance 
 * @param [in] tileopcounters:     count the number of tiles 
 * @param [in] opcounters:         count operations during factorization 
 * @param [in] critical_path_time: measure the critical path time 
 * @param [intout] info            0 on all nodes if successful. 
 *                                 > 0 if the leading minor of order i of A is not positive
 *                                 definite, so the factorization could not be completed, and the
 *                                 solution has not been computed. Info will be equal to i on the
 *                                 node that owns the diagonal element (i,i), and 0 on all other nodes
 *
 */
parsec_taskpool_t*
HiCMA_dpotrf_L_2flow_New( parsec_context_t *parsec,
                  int uplo,
                  parsec_tiled_matrix_dc_t *A,
                  parsec_tiled_matrix_dc_t *Ar,
                  parsec_tiled_matrix_dc_t *RG,
                  parsec_tiled_matrix_dc_t *Rank,
                  double acc,
		  int fixed_rank,
                  int maxrank,
		  int lookahead,
                  int band_size,
                  int hmb,
		  int compmaxrank,
                  int send_full_tile,
                  unsigned long* tileopcounters,
                  unsigned long* opcounters,
                  double *critical_path_time,
                  int *info );

/* 2flow version, destructor */
void HiCMA_dpotrf_L_2flow_Destruct(parsec_taskpool_t* _tp);

/* 3flow version non-blocking 
 *
 * @param [in] parsec:             parsec context 
 * @param [inout] A:               the data, already distributed and allocated
 * @param [inout] Au:              the data, already distributed and allocated
 * @param [inout] Av:              the data, already distributed and allocated
 * @param [inout] Ar:              the rank info, already distributed and allocated
 * @param [in] RG:                 used for reordering GEMM 
 * @param [in] Rank:               used for gathering time info during factorization 
 * @param [in] acc:                accuracy threshold 
 * @param [in] fixed_rank:         fixed rank threshold used in recompression stage of HCORE_GEMM 
 * @param [in] maxrank:            max rank threshold for storage
 * @param [in] lookahead:          lookahead to guide workflow 
 * @param [in] band_size:          band size to control dense tiles; band_size=1: only diagonal in dense 
 * @param [in] hmb:                hierarchical mb; the sub tile size used in the recursive kernel 
 * @param [in] compmaxrank:        max rank threshold used in computation 
 * @param [in] send_full_tile:     whether send full tile during factorization; default is false to give better performance 
 * @param [in] tileopcounters:     count the number of tiles 
 * @param [in] opcounters:         count operations during factorization 
 * @param [in] critical_path_time: measure the critical path time 
 * @param [intout] info            0 on all nodes if successful. 
 *                                 > 0 if the leading minor of order i of A is not positive
 *                                 definite, so the factorization could not be completed, and the
 *                                 solution has not been computed. Info will be equal to i on the
 *                                 node that owns the diagonal element (i,i), and 0 on all other nodes
 *
 */
parsec_taskpool_t*
HiCMA_dpotrf_L_3flow_New( parsec_context_t *parsec,
                  int uplo,
                  parsec_tiled_matrix_dc_t *A,
                  parsec_tiled_matrix_dc_t *Au,
                  parsec_tiled_matrix_dc_t *Av,
                  parsec_tiled_matrix_dc_t *Ar,
                  parsec_tiled_matrix_dc_t *RG,
                  parsec_tiled_matrix_dc_t *Rank,
                  double acc,
		  int fixed_rank,
                  int maxrank,
		  int lookahead,
                  int band_size,
                  int hmb,
		  int compmaxrank,
                  int send_full_tile,
                  unsigned long* tileopcounters,
                  unsigned long* opcounters,
                  double *critical_path_time,
                  int *info );

/* 3flow version, destructor */
void HiCMA_dpotrf_L_3flow_Destruct(parsec_taskpool_t* _tp);

#endif /* _HICMA_PARSEC_INTERNAL_H */