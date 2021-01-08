/**
 * @copyright (c) 2020 King Abdullah University of Science and Technology (KAUST).
 *                     The Universiy of Tennessee and The Universiy of Tennessee Research Foundation.
 *                     All rights reserved.
 **/
#include "hicma_parsec.h"

#define FLOAT   double
#define GEMM    dgemm_
void GEMM(char *, char *, int *, int *, int *, FLOAT *, FLOAT *, int *, FLOAT *, int *, FLOAT *, FLOAT *, int *);

#define LOOP    4

int main(int argc, char ** argv)
{
    MPI_Init(NULL, NULL);
    //printf("testing_dense_gemm\n");
    int m, i, j, tlr_rk = -1;
    FLOAT *a, *b, *c;
    FLOAT alpha = 1.;
    FLOAT beta  = 1.;
    double gflops, dstart, dstop;

    if (argc <= 1) {
        printf("Please specify problem size.\n");
        exit(1);
    }

    int type = atoi(argv[1]);
    m = atoi(argv[2]);
    if(argc == 4)
        tlr_rk = atoi(argv[3]);
    else
        tlr_rk = (int)(m/2);
    
    if(type == 0) { /* run dense version */
	    gflops = 2. * (double)m * (double)m * (double)m;

	    a = (FLOAT *)memalign(128, (size_t)m * (size_t)m * sizeof(FLOAT));
	    if (a == (FLOAT *)NULL) {
		    printf("Out of memory for A matrix\n");
		    exit(1);
	    }
	    b = (FLOAT *)memalign(128, (size_t)m * (size_t)m * sizeof(FLOAT));
	    if (b == (FLOAT *)NULL) {
		    printf("Out of memory for B matrix\n");
		    exit(1);
	    }
	    c = (FLOAT *)memalign(128, (size_t)m * (size_t)m * sizeof(FLOAT));
	    if (c == (FLOAT *)NULL) {
		    printf("Out of memory for C matrix\n");
		    exit(1);
	    }

	    for (j = 0; j < m; j ++) {
		    for (i = 0; i < m; i ++) {
			    a[i + (size_t)j * (size_t)m] = rand() / RAND_MAX - 0.5;
			    b[i + (size_t)j * (size_t)m] = rand() / RAND_MAX - 0.5;
			    c[i + (size_t)j * (size_t)m] = rand() / RAND_MAX - 0.5;
		    }
	    }

	    GEMM("N", "N", &m, &m, &m, &alpha, a, &m, b, &m, &beta, c, &m);

	    dstart = MPI_Wtime();

	    for (i = 0; i < LOOP; i ++) {
		    GEMM("N", "N", &m, &m, &m, &alpha, a, &m, b, &m, &beta, c, &m);
	    }

	    dstop = MPI_Wtime();
	    printf("DGEMM Performance N = %6d : %10.4f GF, time %10.6f (s)\n", m, gflops / (dstop - dstart) * (double)LOOP * 1.e-9, (dstop - dstart) / (double)LOOP);
    } else {
	    /* Benchmarking TLR GEMM  */
	    //printf("testing_tlr_gemm\n");
	    FLOAT *Au, *Bu, *Cu;
	    FLOAT *Av, *Bv, *Cv;
	    int *Ar, *Br, *Cr;

	    Au = (FLOAT *)memalign(128, (size_t)m * (size_t)m * sizeof(FLOAT));
	    Av = (FLOAT *)memalign(128, (size_t)m * (size_t)m * sizeof(FLOAT));
	    Bu = (FLOAT *)memalign(128, (size_t)m * (size_t)m * sizeof(FLOAT));
	    Bv = (FLOAT *)memalign(128, (size_t)m * (size_t)m * sizeof(FLOAT));
	    Cu = (FLOAT *)memalign(128, (size_t)m * (size_t)m * sizeof(FLOAT));
	    Cv = (FLOAT *)memalign(128, (size_t)m * (size_t)m * sizeof(FLOAT));
	    Ar = malloc(sizeof(int));
	    Br = malloc(sizeof(int));
	    Cr = malloc(sizeof(int));
	    int ldamu = m;
	    int ldamv = m;
	    int rk = 0;
	    int storage_maxrank = (int)(m/2);
	    int computation_maxrank = (int)(m/2);
	    double acc = 1.0e-8;
	    double* p_elem_work = malloc(m*m*sizeof(double));
        int info = 0;

	    /* init values to the specified tlk_rk based on TMS routine */
	    *Ar = tlr_rk;
	    *Br = tlr_rk;
	    *Cr = tlr_rk;

        /** Temporary U and V for storing output of SVD and transpose V */
        size_t nelm_tUV = m * m;
        int nelm_d = m;
        double* _tU = calloc(nelm_tUV, sizeof(*_tU));
        double* _tV = calloc(nelm_tUV, sizeof(*_tV));
        double* _tUV = calloc(nelm_tUV, sizeof(*_tUV));
        int iseed[4] = {0, 0, 0, 1};
        double * matrix = malloc(m*m*sizeof(double));
        double* _d = calloc(m, sizeof(*_d));
        for(int i = 0; i < m; i++){
            if(i < tlr_rk)
                _d[i] = pow(10, -1*3);
            else
                _d[i] = pow(10, -1*15);
        }
        double* _S = calloc(nelm_d, sizeof(*_d)); 
        double* superb = calloc(nelm_d, sizeof(*_d)); 
        info = LAPACKE_dlatms(LAPACK_COL_MAJOR, m, m, 'U', iseed, 'N', _d, 0, -1, -1, m-1, m-1, 'N', matrix, m);
        if(info != 0) {printf("Error in LAPACKE_dlatms. Info=%d\n", info);return -1;}
        /** Compress */
        /** https://software.intel.com/en-us/node/521150 */
        info = LAPACKE_dgesvd(LAPACK_COL_MAJOR, 'A', 'A', m, m, matrix, m, _S, _tU, m, _tV, m, superb);
        if(info != 0) {printf("Error in LAPACKE_dgesvd. Info=%d\n", info);return -1;}
        /** S * V */
        for(int k = 0; k < m /** ? */; k++){
            double diagval = _S[k];
            cblas_dscal(m, diagval, &_tV[k], m);
        }
        char chall = 'A';
        dlacpy_(&chall, &m, &tlr_rk, _tU, &m, Au, &m);
        LAPACKE_dge_trans(LAPACK_COL_MAJOR, tlr_rk, m, _tV, m, Av, m);

        dlacpy_(&chall, &m, &tlr_rk, _tU, &m, Bu, &m);
        LAPACKE_dge_trans(LAPACK_COL_MAJOR, tlr_rk, m, _tV, m, Bv, m);
        
        dlacpy_(&chall, &m, &tlr_rk, _tU, &m, Cu, &m);
        LAPACKE_dge_trans(LAPACK_COL_MAJOR, tlr_rk, m, _tV, m, Cv, m);

	    /* warm up run */
	    HCORE_dgemm( PlasmaNoTrans, PlasmaTrans, // TODO this call is not ready yet
			    m, /* mb */ 
			    m, /* nb */
			    (double)-1.0, 
			    Au, Av, Ar, ldamu,
			    Bu, Bv, Br, ldamv,
			    (double)1.0,
			    Cu, Cv, Cr, ldamu,
			    rk, storage_maxrank, computation_maxrank, acc, p_elem_work);

	    dstart = MPI_Wtime();
	    for (i = 0; i < LOOP; i ++) {
		    *Ar = tlr_rk;
		    *Br = tlr_rk;
		    *Cr = tlr_rk;
		    HCORE_dgemm( PlasmaNoTrans, PlasmaTrans, // TODO this call is not ready yet
				    m, /* mb */ 
				    m, /* nb */
				    (double)-1.0, 
				    Au, Av, Ar, ldamu,
				    Bu, Bv, Br, ldamv,
				    (double)1.0,
				    Cu, Cv, Cr, ldamu,
				    rk, storage_maxrank, computation_maxrank, acc, p_elem_work);
	    }
	    dstop = MPI_Wtime();

	    printf("TLR DGEMM Performance N = %6d maxrank = %d %10.4f GF, time %10.6f (s)\n", m, tlr_rk, (36. * m * tlr_rk * tlr_rk + 157 * tlr_rk * tlr_rk) / (dstop - dstart) * (double)LOOP * 1.e-9, (dstop - dstart) / (double)LOOP);
    }
    MPI_Finalize();
}
