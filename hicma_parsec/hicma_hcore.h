/**
 * @copyright (c) 2020 King Abdullah University of Science and Technology (KAUST).
 *                     The Universiy of Tennessee and The Universiy of Tennessee Research Foundation.
 *                     All rights reserved.
 **/
#ifndef _HICMA_HCORE_H_
#define _HICMA_HCORE_H_
//#include "plasma.h"
#define HiCMA_HCORE_ERROR_THRESHOLD_DENSE 1e-14
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
void HiCMA_HCORE_dgemm(int transA, int transB,
        int M, int N,
        double alpha,
        double *AU,
        double *AV,
        int *Ark,
        int LDA,
        double *BU,
        double *BV,
        int *Brk,
        int LDB,
        double beta,
        double *CU,
        double *CV,
        int *Crk,
        int LDC,
        int rk,
        int storage_maxrk,
        int maxrk,
        double acc,
        double* work
);
void HiCMA_HCORE_dgemm_qr_svd_b_dense(int transA, int transB,
        int M, int N,
        double alpha,
        double *AU,
        double *AV,
        int *Ark,
        int LDA,
        double *B,
        int LDB,
        double beta,
        double *CU,
        double *CV,
        int *Crk,
        int LDC,
        int rk,
        int storage_maxrank,
        int maxrk, /*this is compmaxrank*/
        double acc,
        double* work,
        /** parameters that will be passed to HiCMA_HCORE_dgemm_ormqr */
        double*      *_p_work_new     ,
        double*      *_p__CU          ,
        double*      *_p__CV          ,
        int          *_p_CU_ncols     ,
        int          *_p_new_UVrk     ,
        double*      *_p_newU         ,
        int          *_p_ld_newU      ,
        double*      *_p_qrtauA       ,
        int          *_p_CV_ncols     ,
        double*      *_p_newV         ,
        int          *_p_ld_newV      ,
        double*      *_p_qrtauB       ,
        int          *_p_use_CUV_clone,
        double*      *_p_CUclone      ,
        int          *_p_ld_CUclone   ,
        double*      *_p__CU_save     ,
        double*      *_p_CVclone      ,
        int          *_p_ld_CVclone   ,
        double*      *_p__CV_save     
);
void HiCMA_HCORE_dgemm_qr_svd(int transA, int transB,
        int M, int N,
        double alpha,
        double *AU,
        double *AV,
        int *Ark,
        int LDA,
        double *BU,
        double *BV,
        int *Brk,
        int LDB,
        double beta,
        double *CU,
        double *CV,
        int *Crk,
        int LDC,
        int rk,
        int storage_maxrank,
        int maxrk, /*this is compmaxrank*/
        double acc,
        double* work,
        /** parameters that will be passed to HiCMA_HCORE_dgemm_ormqr */
        double*      *_p_work_new     ,
        double*      *_p__CU          ,
        double*      *_p__CV          ,
        int          *_p_CU_ncols     ,
        int          *_p_new_UVrk     ,
        double*      *_p_newU         ,
        int          *_p_ld_newU      ,
        double*      *_p_qrtauA       ,
        int          *_p_CV_ncols     ,
        double*      *_p_newV         ,
        int          *_p_ld_newV      ,
        double*      *_p_qrtauB       ,
        int          *_p_use_CUV_clone,
        double*      *_p_CUclone      ,
        int          *_p_ld_CUclone   ,
        double*      *_p__CU_save     ,
        double*      *_p_CVclone      ,
        int          *_p_ld_CVclone   ,
        double*      *_p__CV_save     
);
void HiCMA_HCORE_dgemm_ormqr(int transA, int transB,
        int M, int N,
        double alpha,
        double *AU,
        double *AV,
        int *Ark,
        int LDA,
        double beta,
        double *CU,
        double *CV,
        int *Crk,
        int LDC,
        int rk,
        int storage_maxrank,
        int maxrk, /*this is compmaxrank*/
        double acc,
        double* work,
        /** parameters coming from HiCMA_HCORE_dgemm_qr_svd */
        double* _CU,
        double* _CV,
        int CU_ncols,
        int new_UVrk,
        double* newU,
        int ld_newU,
        double* qrtauA,
        int CV_ncols,
        double* newV,
        int ld_newV,
        double* qrtauB,
        int use_CUV_clone,
        double* CUclone,
        int ld_CUclone,
        double* _CU_save,
        double* CVclone,
        int ld_CVclone,
        double* _CV_save
);
/*
 * Copies _B_ into p_elem_work and returns InfNorm(B-A)
 */
double
HiCMA_check_norm(int _n_, int _ld_, double* _A_, double* _B_, double* p_elem_work);
#endif
