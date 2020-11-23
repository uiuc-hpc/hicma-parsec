/**
 * @copyright (c) 2020 King Abdullah University of Science and Technology (KAUST).
 *                     The Universiy of Tennessee and The Universiy of Tennessee Research Foundation.
 *                     All rights reserved.
 **/
#ifndef _HICMA_PARSEC_H_
#define _HICMA_PARSEC_H_

#include "hicma_parsec_internal.h"

/* Main routine for TLR Cholesky
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
 * @param [in] two_flow:           force to run two_flow version even if band_size == 1
 * @param [in] tileopcounters:     count the number of tiles 
 * @param [in] opcounters:         count operations during factorization 
 * @param [in] critical_path_time: measure the critical path time 
 * @return info:                   0 on all nodes if successful. 
 *                                 > 0 if the leading minor of order i of A is not positive
 *                                 definite, so the factorization could not be completed, and the
 *                                 solution has not been computed. Info will be equal to i on the
 *                                 node that owns the diagonal element (i,i), and 0 on all other nodes
 */
int HiCMA_dpotrf_L( parsec_context_t* parsec,
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
              int two_flow, 
              unsigned long* tileopcounters,
              unsigned long* opcounters,
              double *critical_path_time
              );

/**
 * Set offset in Av based on A for the 3flow version
 *
 * @param [in] dcA:          the data, already distributed and allocated
 * @param [inout] dcAv:      the data, already distributed and allocated
 * @param [in] dcAr:         the rank data, already distributed and allocated
 * @param [in] maxrank:      max rank
 */
int parsec_Av_memory(parsec_context_t *parsec,
                     parsec_tiled_matrix_dc_t *dcA,
                     parsec_tiled_matrix_dc_t *dcAv,
                     parsec_tiled_matrix_dc_t *dcAr,
                     int maxrank);

/**
 * @brief free memory
 *
 * @param [inout] dcA: the data, already distributed and allocated
 * @param [in] band_size:   band size
 * @param [in] indicator:   0, only band; otherwise all
 */
int parsec_band_free(parsec_context_t *parsec,
                     parsec_tiled_matrix_dc_t *dcA,
                     int band_size, int indicator);

/**
 * Generate tile on band
 *              
 * @param [in] parsec:       parsec context
 * @param [in] uplo:         support PlasmaLower now
 * @param [inout] dcA:       the data, already distributed and allocated
 * @param [in] data:         data in generate tiles
 * @param [in] kernel:       kernels in generate matrix
 * @param [in] index:        index in generate matrix
 * @param [in] band_size:    band size
 */
int parsec_band_regenerate( parsec_context_t *parsec,
                int uplo,
                parsec_tiled_matrix_dc_t *dcA,
                void *data,
                STARSH_kernel *kernel,
                STARSH_int *index,
                int band_size);

/**
 * Band size auto tuning
 *
 * @param [in] parsec:       parsec context
 * @param [in] dcAr:         the data, already distributed and allocated
 * @param [in] dcFake:       used for distribution
 * @param [in] Ar_copy:      Rank info
 * @param [in] NB:           tile size
 * @param [in] disp:         displacement for each process in memory
 * @param [in] nb_elem_r:    number of elements in row
 */
int parsec_band_size_auto_tuning(parsec_context_t *parsec,
                          parsec_tiled_matrix_dc_t *dcAr,
                          parsec_tiled_matrix_dc_t *dcFake,
                          int *Ar_copy,
                          int NB,
                          int *disp,
                          int *nb_elem_r);

/**
 * @brief Gather dcY to rank_array
 *
 * @param [in] dcY: the data, already distributed and allocated
 * @param [in] rank_array: array of rank
 */
int parsec_rank_gather(parsec_context_t *parsec,
                parsec_tiled_matrix_dc_t *dcY,
                int *rank_array);

/**
 * Check rank correctness and set -1 to tiles on band
 *
 * @param [in] parsec:       parsec context 
 * @param [inout] dcAr:      the rank data, already distributed and allocated
 * @param [in] band_size:    band size
 */
int parsec_rank_check(parsec_context_t *parsec,
                      parsec_tiled_matrix_dc_t *dcAr,
                      int band_size);

/**
 * Gather and print rank distribution
 *
 * @param [in] parsec:       parsec context
 * @param [inout] dcRank:    the rank data, already distributed and allocated
 * @param [in] band_size:    band size
 */
int parsec_rank_print(parsec_context_t *parsec,
                      parsec_tiled_matrix_dc_t *dcRank,
                      int band_size);

/** Uncompresses approximate matrix dcA into dcA0 
 *
 * @param [in] parsec:       parsec context
 * @param [in] uplo:         support PlasmaLower now
 * @param [inout] dcA0:      the data, already distributed and allocated
 * @param [in] dcA:          the data, already distributed and allocated
 * @param [in] dcAr:         the rank data, already distributed and allocated
 * @param [in] band_size:    band size
 * @param [inout] info:      info to check result
 */
int STARSH_check( parsec_context_t *parsec,
                int uplo,
                parsec_tiled_matrix_dc_t *dcA0,
                parsec_tiled_matrix_dc_t *dcA,
                parsec_tiled_matrix_dc_t *dcAr,
                int band_size,
                int *info);

/**
 * Generate matrix
 *
 * @param [in] parsec:       parsec context
 * @param [in] uplo:         support PlasmaLower now
 * @param [inout] dcA:       the data, already distributed and allocated
 * @param [in] dcAr:         the rank data, already distributed and allocated
 * @param [in] data:         data in generate tiles
 * @param [in] kernel:       kernels in generate matrix
 * @param [in] index:        index in generate matrix
 * @param [in] tol:          fixed accuracy threshold
 * @param [in] maxrank:      max rank
 * @param [in] band_size:    band size
 * @param [in] send_full_tile: when send full tile
 * @param [inout] info:        check result
 */
int STARSH_gen( parsec_context_t *parsec,
                int uplo,
                parsec_tiled_matrix_dc_t *dcA,
                parsec_tiled_matrix_dc_t *dcAr,
                void *data,
                STARSH_kernel *kernel,
                STARSH_int *index,
                double tol,
                int maxrank,
                int band_size,
                int send_full_tile,
                int *info);

/* Set up parsec*/
parsec_context_t *setup_parsec(int argc, char* argv[], int *iparam, double *dparam);

/* Clean parsec */
void cleanup_parsec(parsec_context_t* parsec, int *iparam);

#if defined(PARSEC_HAVE_CUDA)
typedef struct parsec_potrf_stream_workspace_s {
    cusolverDnHandle_t handle;
    void *gpu_buffer;
    int buffer_size;
}parsec_potrf_stream_workspace_t;

typedef struct parsec_potrf_gpu_workspace_s {
    parsec_potrf_stream_workspace_t *stream_workspace; 
    parsec_device_cuda_module_t *gpu_device; 
}parsec_potrf_gpu_workspace_t;

typedef struct parsec_potrf_workspace_s {
    parsec_potrf_gpu_workspace_t *gpu_workspace;
    int info; 
}parsec_potrf_workspace_t;
#endif

#endif
