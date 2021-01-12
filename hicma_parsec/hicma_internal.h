/**
 * @copyright (c) 2020 King Abdullah University of Science and Technology (KAUST).
 *                     The Universiy of Tennessee and The Universiy of Tennessee Research Foundation.
 *                     All rights reserved.
 **/
#ifndef _HICMA_INTERNAL_H_
#define _HICMA_INTERNAL_H_
//#include "plasma.h"
#include "hicma_struct.h"
#include "hicma_parsec_internal.h"

int HICMA_init();
int HICMA_get_print_mat();
int HICMA_set_print_mat();
int HICMA_get_print_index();
int HICMA_get_print_index_end();
int HICMA_set_print_index();
int HICMA_set_print_index_end();
int HICMA_unset_print_index_end();
int HICMA_set_use_fast_hcore_zgemm();
int HICMA_set_starsh_format(STARSH_blrf *starsh_format);
STARSH_blrf* HICMA_get_starsh_format();
int HICMA_get_use_fast_hcore_zgemm();
int HICMA_get_always_fixed_rank();
int HICMA_get_fixed_rank();
int HICMA_set_fixed_rank(int rank);
void HICMA_get_stat(char uplo, int *Ark, size_t m, size_t n, size_t ld,  HICMA_stat_t *stat);
void HICMA_get_stat2(int *Ark, size_t m, int maxrank,  HICMA_stat_t *stat);
void HICMA_print_stat(HICMA_stat_t stat);
int
tile_dpotrf( int /*int*/ uplo,
                   int m, 
                   double* A, int lda
                   , int Am, int An // tile index, might be removed
                   ) ;
int
tile_dtrsm( int /*int*/ side,
                   int /*int*/ uplo,
                   int /*int*/ transa,
                   int /*int*/ diag,
                   int m, int n,
                   double alpha,
                   double* A, int lda,
                   double* B, int ldb,
                   int* Brk, int Am, int An, int Bm, int Bn
                   );
#endif
