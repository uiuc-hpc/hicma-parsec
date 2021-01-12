/**
 * @copyright (c) 2020 King Abdullah University of Science and Technology (KAUST).
 *                     The Universiy of Tennessee and The Universiy of Tennessee Research Foundation.
 *                     All rights reserved.
 **/
#include "hicma_internal.h"
int use_scratch = 1;
/*
mpicc `pkg-config --cflags plasma` `pkg-config --libs plasma` -c hicma_hcore.c  && ar rc libhicmahcore.a hicma_hcore.o
 */
int trsm_print_index_end = 0;
int
tile_dpotrf( int /*int*/ uplo,
                   int m, 
                   double* A, int lda, int Am, int An) {
    /*printf("===========%s %s %d==========\n", __FILE__, __func__, __LINE__);*/
    if(HICMA_get_print_index() == 1){
        printf("%d+POTRF\t|AD(%d,%d) m:%d lda(11):%d\n", 0, Am, An, m, lda);
    }
    if(HICMA_get_print_mat() == 1){
        printf("%d\tpotrf-input A\n", __LINE__);
        _printmat(A, m, m, lda);
    }
    int iinfo = 0;

#if 0
        printf("\n before dpotrf(%d)\n", Am);
        for(int i = 0; i < 2; i++) { 
            for(int j = 0; j < 2; j++) { 
                printf("%lf ", A[j*lda+i]);
            }
            printf("\n");
         }
        printf("\n");
#endif


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
int
tile_dtrsm( int /*int*/ side,
                   int /*int*/ uplo,
                   int /*int*/ transA,
                   int /*int*/ diag,
                   int m, int n,
                   double alpha,
                   double* A, int lda,
                   double* B, int ldb,
                   int* Brk, int Am, int An, int Bm, int Bn
                   ) {
    /*printf("%s %s %d\n", __FILE__, __func__, __LINE__);*/
    //return 0;
    struct timeval tvalBefore, tvalAfter;  // removed comma
    gettimeofday (&tvalBefore, NULL);
    //int m;
    double *BUV;
    /*int nBUV;*/

    //A = (double *)STARPU_MATRIX_GET_PTR(descr[0]);
    //BUV = (double *)STARPU_MATRIX_GET_PTR(descr[1]);
    //Brk = (double *)STARPU_MATRIX_GET_PTR(descr[2]);
    int _Brk;
    _Brk = Brk[0];

    /*int nBU = nBUV/2;*/
    /*size_t nelm_BU = (size_t)ldb * (size_t)nBU;*/
    /*double *B = &(BUV[nelm_BU]);*/

    /*CORE_ztrsm(side, uplo,*/
        /*transA, diag,*/
        /*m, n,*/
        /*alpha, A, lda,*/
        /*B, ldb);*/
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
    if(HICMA_get_print_index() == 1 || HICMA_get_print_index_end() == 1 || trsm_print_index_end){
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


void HICMA_get_stat(char uplo, int *Ark, size_t m, size_t n, size_t ld,  HICMA_stat_t *stat)
{
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

            if(0){
                if(jnt<imt)
                    printf("Tile %d,%d %d\n", imt, jnt, rank);
            }
        }
    }
    final_avgrank = (final_totalrank)/(ntiles*1.0);
    stat->min = minrank;
    stat->max = final_maxrank;
    stat->avg = final_avgrank;
}

/* Used for Gather rank of low-rank tiles */ 
void HICMA_get_stat2(int *G, size_t lda, int band_size, HICMA_stat_t *stat)
{
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

void HICMA_print_stat(HICMA_stat_t stat)
{
    printf("avg:%g min:%d max:%d\n", stat.avg, stat.min, stat.max);
}
