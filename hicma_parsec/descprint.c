/**
 * @copyright (c) 2021     King Abdullah University of Science and Technology (KAUST).
 * @copyright (c) 2021     The Universiy of Tennessee and The Universiy of Tennessee Research Foundation.
 *                         All rights reserved.
**/
/**
 * @file descprint.c
 *
 * This file contains the functions for printing matrices 
 * 
 * @version 0.1.0
 * @date 2021-01-24
 **/
#include "auxdescutil.h"
int _nelm_limit = 1200;
int _nrows_limit = 37;
int _ncols_limit = 33;
int _ndim_limit = 37;
void _printmat(double * A, int m, int n, int ld){
    printf("%s@%d M:%d N:%d LD:%d %p [\n", __FILE__, __LINE__, m, n, ld, A);
    int i, j, nelm = 0;
    for(i=0;i<m;i++){
        printf("[");
        for(j=0;j<n;j++){
            printf("%+.4e", A[j*ld+i]);
            if(j!=n-1){
                printf(",");
            }
            nelm++;
            if(nelm >= _nelm_limit){
                printf("\n");
                return;
            }
            if(j==_ncols_limit)
                break;
        }
        printf("]");
        if(i!=m-1){
            printf(",");
            printf("\n");
        }
        if(i==_nrows_limit)
            break;
    }
    printf("]\n");
}

