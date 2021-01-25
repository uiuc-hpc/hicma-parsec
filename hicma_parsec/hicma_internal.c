/**
 * @copyright (c) 2021 King Abdullah University of Science and Technology (KAUST).
 *                     The Universiy of Tennessee and The Universiy of Tennessee Research Foundation.
 *                     All rights reserved.
 *
 * This file declares the HICMA auxiliary functions. 
 *
 * @version 0.1.0
 * @date 2021-01-24
 **/
#include "hicma_internal.h"

// Use work array while using HCORE functions.
int use_scratch = 1;

int tile_dpotrf( int uplo, int m, double* A, int lda, int Am, int An) {
    if(HICMA_get_print_index() == 1){
        printf("%d+POTRF\t|AD(%d,%d) m:%d lda(11):%d\n", 0, Am, An, m, lda);
    }
    if(HICMA_get_print_mat() == 1){
        printf("%d\tpotrf-input A\n", __LINE__);
        _printmat(A, m, m, lda);
    }
    int iinfo = 0;

    CORE_dpotrf( uplo, m, A, lda, &iinfo );

    if(iinfo != 0){
        printf("%s %d: dpotrf failed with a return value of %d. uplo:%d m:%d A:%p lda:%d\n", __FILE__, __LINE__, iinfo, uplo, m, A, lda);
        fflush(stdout);
        exit(-1);
    }
    if(HICMA_get_print_mat() == 1){
        printf("%d\tpotrf-output A\n", __LINE__);
        _printmat(A, m, m, lda);
    }
    return iinfo;
}
int tile_dtrsm( int side, int uplo, int transA, int diag, int m, int n, double alpha, double* A, int lda, double* B, int ldb, int* Brk, int Am, int An, int Bm, int Bn) {
    struct timeval tvalBefore, tvalAfter;  // removed comma
    gettimeofday (&tvalBefore, NULL);
    double *BUV;
    int _Brk;
    _Brk = Brk[0];
    if(HICMA_get_print_index() == 1){
        printf("%d+TRSM\t|AD(%d,%d) BV(%d,%d)%d m:%d lda(11):%d ldb(12):%d\n", 0,Am,An, Bm, Bn, _Brk, m, lda, ldb);
    }
    if(HICMA_get_print_mat() == 1){
        printf("%d\ttrsm-input A\n", __LINE__);
        _printmat(A, m, m, lda);
        printf("%d\ttrsm-input B\n", __LINE__);
        _printmat(B, m, _Brk, ldb);
    }
    cblas_dtrsm(
        CblasColMajor,
        side, uplo,
        transA, diag,
        m,
        _Brk,
        alpha, A, lda,
        B, ldb);
    if(HICMA_get_print_index() == 1 || HICMA_get_print_index_end() == 1){
        gettimeofday (&tvalAfter, NULL);
        printf("%d-TRSM\t|AD(%d,%d)%dx%d-%d BV(%d,%d)%dx%d-%d m:%d\t\t\t\tTRSM: %.4f\n", 0, Am,An, m, m, lda,Bm, Bn, m, _Brk, ldb, m,
                (tvalAfter.tv_sec - tvalBefore.tv_sec)
                 +(tvalAfter.tv_usec - tvalBefore.tv_usec)/1000000.0
                );
    }
    if(HICMA_get_print_mat() == 1){
        printf("%d\ttrsm-output\n", __LINE__);
        _printmat(B, m, _Brk, ldb);
    }
    return 0;
}

void HICMA_get_stat(char uplo, int *Ark, size_t m, size_t n, size_t ld,  HICMA_stat_t *stat) {
    double final_avgrank;
    int final_maxrank = 0;
    int minrank = 10000;
    int final_totalrank = 0;
    int *MAT = Ark;
    int64_t i, j, imt, jnt, nelm = 0;
    int ntiles = 0;

    for(imt=0;imt<m;imt++){
        for(jnt=0;jnt<n;jnt++){
            if(imt == jnt)
                continue;
            if(uplo == 'L' && imt < jnt)
                continue;
            if(uplo == 'U' && imt > jnt)
                continue;
            int *A = MAT+imt+jnt*ld;
            int rank = A[0];
            if(rank > final_maxrank){
                final_maxrank = rank;
            }
            if(rank < minrank){
                minrank = rank;
            }
            final_totalrank += rank;
            ntiles++;
        }
    }
    final_avgrank = (final_totalrank)/(ntiles*1.0);
    stat->min = minrank;
    stat->max = final_maxrank;
    stat->avg = final_avgrank;
}

/* Used for Gather rank of low-rank tiles */ 
void HICMA_get_stat2(int *G, size_t lda, int band_size, HICMA_stat_t *stat) {
    int min = INT_MAX, max = 0;
    long int num = 0;
    long long int sum = 0;

    for(int i = band_size; i < lda; i++){
        for(int j = 0; j < i-band_size+1; j++){
            sum += G[j*lda+i];
	    num++;
            if( G[j*lda+i] < min )
                min = G[j*lda+i];
            if( G[j*lda+i] > max )
                max = G[j*lda+i];
        } 
    }

    stat->min = min;
    stat->max = max;
    stat->avg = ((double)sum)/num;
}

void HICMA_print_stat(HICMA_stat_t stat) {
    printf("avg:%g min:%d max:%d\n", stat.avg, stat.min, stat.max);
}
