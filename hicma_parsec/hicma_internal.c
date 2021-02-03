/**
 * @copyright (c) 2021     King Abdullah University of Science and Technology (KAUST).
 * @copyright (c) 2021     The Universiy of Tennessee and The Universiy of Tennessee Research Foundation.
 *                         All rights reserved.
 *
 * This file declares the HICMA auxiliary functions. 
 *
 * @version 0.1.0
 * @date 2021-01-24
 **/
#include "hicma_internal.h"

// Use work array while using HCORE functions.
int use_scratch = 1;

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
