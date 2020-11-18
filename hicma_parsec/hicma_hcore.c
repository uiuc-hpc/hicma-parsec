/**
 * @copyright (c) 2020 King Abdullah University of Science and Technology (KAUST).
 *                     The Universiy of Tennessee and The Universiy of Tennessee Research Foundation.
 *                     All rights reserved.
 **/
#include "hicma_hcore.h"

/*
mpicc `pkg-config --cflags plasma` `pkg-config --libs plasma` -c hicma_hcore.c  && ar rc libhicmahcore.a hicma_hcore.o
 */
int trsm_print_index_end = 0;
int
HiCMA_HCORE_dpotrf( int /*int*/ uplo,
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
HiCMA_HCORE_dtrsm( int /*int*/ side,
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
int syrk_print_index = 0;
int use_scratch = 1;

/*
 * N != LDA
 * Leading dim of AU is LDA
 *                                                            M LDA
 *                                                          ______
 *                                                         |  AU  |
 *                                                     K   |______|
 *                                                    ____    
 *                                                   |    |            
 *                                                 N | AV |
 *                                               LDA |    |
 *                                                   |____|            
 *                                       N LDA    X
 *                                      _________
 *                                     |   AV    |
 *         M                        K  |_________|
 *       _____      _____         ____    
 *      |     |    |     |     M |    |               
 *    M | CD  | =  | CD  | + LDA | AU |             
 * LDCD |_____|    |_____|       |____|             
 */
void HiCMA_HCORE_dsyrk(int uplo, int trans,
                int M, int K,
                double alpha,
                const double *AU, int LDAU,
                const double *AV, int LDAV,
                double beta, 
                double *CD, int LDCD,
                double* work
                )
{
    int64_t N = LDAU; //ASSUMPTION FIXME
    int64_t LDA = LDAU;

    /*cblas_zsyrk(*/
            /*CblasColMajor,*/
            /*(CBLAS_UPLO)uplo, (CBLAS_TRANSPOSE)trans,*/
            /*N, N,*/
            /*CBLAS_SADDR(alpha), A, LDA,*/
            /*CBLAS_SADDR(beta), C, LDC);*/

    int64_t A_colfactor_ncols = K;
    /// C = C + alpha * A * A'
    /// C = C + alpha * ( (A^u * (A^v * A^v^T) ) * A^u^T)
    /// A^v * B^v^T
    int64_t bufmtx_nrows = A_colfactor_ncols;
    int64_t bufmtx_ncols = A_colfactor_ncols;
    size_t bufmtx_nelm = bufmtx_nrows * bufmtx_ncols;
    /*size_t bufmtx_nelm = K * K;*/
    //ALLOCATE1 bufmtx_syrk_kk[myid];
    double* bufmtx = NULL;
    if(use_scratch){
        bufmtx = work;
    } else {
        bufmtx = malloc(bufmtx_nelm * sizeof(double));
    }
    assert(bufmtx != NULL);
    double alpha2 = 1.0f;
    double beta2 = 0.0f;
    int64_t ld_AV = M; //ASSUMPTION
    int64_t ld_bufmtx = bufmtx_nrows; //ASSUMPTION

    int64_t AV_nrows = M; //ASSUMPTION
    int64_t AV_ncols = K;
    if(syrk_print_index){
        printf(" SYRK 1 |%d\t|Trans NoTrans AV_ncols:%d AV_ncols:%d AV_nrows:%d alpha2:%.2e ld_AV:%d beta2:%.2e ld_bufmtx:%d AV:%p bufmtx:%p\n", __LINE__,
                AV_ncols, AV_ncols, AV_nrows, alpha2, ld_AV, beta2, ld_bufmtx, AV, bufmtx    
              );
        printf(" SYRK 1 |%d\t|Trans NoTrans K:%d K:%d N:%d alpha2:%.2e LDA:%d beta2:%.2e K:%d AV:%p bufmtx:%p\n", __LINE__,
                                              K, K, N,     alpha2,     LDA,   beta2,     K,   AV, bufmtx
              );
    } 
    /*assert(AV_ncols == K);*/
    /*assert(AV_nrows == N);*/
    /*assert(ld_AV    == LDA);*/
    /*assert(ld_bufmtx== K);*/
    cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, K, K, N, alpha2, AV, LDA, AV, LDA, beta2, bufmtx, K);
    /*bufmtx_nrows = bufmtx_ncols = AV_ncols;  */
    double* AcolBcolT = bufmtx;
    int64_t AcolBcolT_ncols = bufmtx_ncols;

    /// A^u * (A^v * B^v^T)
    int64_t A_rowfactor_nrows = M;
    int64_t A_rowfactor_ncols = K;
    int64_t AU_nrows = M;
    int64_t AU_ncols = K;
    int64_t bufmtx2_nrows = A_rowfactor_nrows; 
    int64_t bufmtx2_ncols = A_rowfactor_ncols;
    //ALLOCATE2
    //bufmtx_syrk_U[myid];
    double* bufmtx2 = NULL;
    
    if(use_scratch){
        bufmtx2 = work + bufmtx_nelm;
    } else {
        /*malloc(bufmtx2_nrows * bufmtx2_ncols * sizeof(double));*/
        bufmtx2 = malloc(M * K * sizeof(double));
    }
    assert(bufmtx2 != NULL);

    int64_t ld_AU = M; //ASSUMPTION
    int64_t ld_bufmtx2 = M; //ASSUMPTION
    if(syrk_print_index){
        printf(" SYRK 2 |%d\t|NoTrans NoTrans AU_nrows:%d AcolBcolT_ncols:%d AU_ncols:%d alpha2:%.2e ld_AU:%d ld_bufmtx:%d beta2:%.2e ld_bufmtx2:%d AU:%p AcolBcolT:%p bufmtx2:%p\n", __LINE__,
                 AU_nrows, AcolBcolT_ncols, AU_ncols, alpha2, ld_AU, ld_bufmtx, beta2, ld_bufmtx2, AU, AcolBcolT, bufmtx2
              );
        printf(" SYRK 2 |%d\t|NoTrans NoTrans M:%d K:%d K:%d alpha2:%.2e LDA:%d K:%d beta2:%.2e M:%d AU:%p AcolBcolT:%p bufmtx2:%p\n", __LINE__,
                                                M,   K,   K,        alpha2, LDA, K,   beta2,     M , AU, AcolBcolT, bufmtx2
              );
        
    } 
    /*assert(AU_nrows           == M);*/
    /*assert(AcolBcolT_ncols    == K);*/
    /*assert(AU_ncols           == K);*/
    /*assert(ld_AU              == LDA);*/
    /*assert(ld_bufmtx2         == M);*/
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, M, K, K, alpha2, AU, LDA, AcolBcolT, K, beta2, bufmtx2, M);
    bufmtx2_nrows = AU_nrows;
    bufmtx2_ncols = AcolBcolT_ncols;
  

    double* Arow_AcolBcolT = bufmtx2; 
    int64_t ld_C = M; //ASSUMPTION
    /// (A^u * (A^v * B^v^T) ) * B^u^T   
    int64_t Arow_AcolBcolT_nrows = bufmtx2_nrows;
    int64_t Arow_AcolBcolT_ncols = bufmtx2_ncols;
    if(syrk_print_index){
        printf(" SYRK 3 |%d\t|NoTrans Trans Arow_AcolBcolT_nrows:%d AU_nrows:%d Arow_AcolBcolT_ncols:%d alpha:%.2e ld_bufmtx2:%d ld_AU:%d alpha2:%.2e LDCD:%d Arow_AcolBcolT:%p AU:%p CD:%p\n", 
                __LINE__,Arow_AcolBcolT_nrows,  AU_nrows,  Arow_AcolBcolT_ncols,  alpha,  ld_bufmtx2,  ld_AU,  alpha2,  LDCD,  Arow_AcolBcolT,  AU,  CD 
              );
        printf(" SYRK 3 |%d\t|NoTrans Trans M:%d M:%d K:%d alpha:%.2e M:%d LDA:%d alpha2:%.2e LDCD:%d Arow_AcolBcolT:%p AU:%p CD:%p\n", 
                __LINE__, M, M,  K,  alpha,  M,  LDA,  alpha2,  LDCD,  Arow_AcolBcolT,  AU,  CD 
              );
    } 
    /*assert(Arow_AcolBcolT_nrows == M);*/
    /*assert(Arow_AcolBcolT_ncols == K);*/
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, 
            M, M, K, 
            alpha, 
            Arow_AcolBcolT, M, 
            AU, LDA, 
            alpha2, 
            CD, LDCD);
    int64_t C_ncols = AU_nrows;
    {
        int i,j;
        for(j = 1; j < M; j++){
            for(i = 0; i < j; i++){
                CD[j*LDCD+i] = beta2;
            }
        }
    }

    if(use_scratch == 0){
        //FREE1
        free(bufmtx);
        //FREE2
        free(bufmtx2);
    }
}

int use_trmm = 1;
int gemm_print_index = 0;
int gemm_print_mat = 0;
int hc_nelm_limit = 512;
void hc_printmat(double * A, int m, int n, int ld){
    printf("M:%d N:%d LD:%d\n[", m, n, ld);
    int i, j, nelm = 0;
    for(i=0;i<m;i++){
        printf("[");
        for(j=0;j<n;j++){
            printf("%+.2e", A[j*ld+i]);
        if(j!=n-1){
            //printf(", ");
        }
            //printf("%g ", A[j*tld(descZ)+i]);
            //printf("%g\t", A[j*descZ->n+i]);
            //printf("(%d,%d,%d) %g\t", i,j,descZ->mb,A[j*descZ->mb+i]);
            //printf("(%d,%d,%d) %g\t", i,j,descZ->n,A[j*descZ->n+i]);
            //printf("(%d,%d,%d) %g [%d %d %d]\t", i,j,descZ->n,A[j*descZ->n+i], descZ->m, descZ->lm, descZ->ln);
            nelm++;
            if(nelm >= hc_nelm_limit){
                printf("\n");
                return;
            }
        }
        printf("]");
        if(i!=m-1){
            printf(",");
            printf("\n");
        }
    }
    printf("]\n");
}
// void __qra(ka_matrix* _CU, ka_matrix* _AU, double alpha, double* qrtauA)
void __qra(int _M,
        int maxrank,
        double* _CU, int ld_CU, int _Crk,
        int* pnew_CU_ncols,
        double* _AU, int ld_AU, int _Ark,
        double alpha, double beta, double* qrtauA){
    int info;
    int AU_nrows = _M;// assume that _M == nrows of A
    int AU_ncols = _Ark;
    int CU_ncols = _Crk;
    if ((_Ark + _Crk) > 2*maxrank){
        fprintf(stderr, "%s %s %d: Sum of ranks (%d) is too big! _Ark:%d _Crk:%d maxrank:%d (x2: %d)\n",
                __FILE__, __func__, __LINE__, (_Ark + _Crk), _Ark, _Crk, maxrank, 2*maxrank);
        exit(-1);
    }
    int nelm_AU = AU_nrows * AU_ncols;
    int incOne = 1;
    //hc_printmat(_CU, _M, _M, ld_CU);
    //hc_printmat(_CU, _M, _M, ld_CU);
    //ERRONEOUS. here is no ld!!!
    if(gemm_print_index){
        printf(" QRA\t|%d\t|nelm_AU:%d alpha:%g CU_ncols:%d ld_CU:%d CU_ncols*ld_CU:%d\n", __LINE__, nelm_AU, alpha, CU_ncols, ld_CU, CU_ncols*ld_CU);
    }
    cblas_dcopy(nelm_AU, _AU, incOne,  &_CU[CU_ncols*ld_CU], incOne);
    //hc_printmat(_CU, _M, _M, ld_CU);
    double d_one = (double)1.0;
    if(alpha != d_one){
        cblas_dscal(nelm_AU, (alpha) , &_CU[CU_ncols*ld_CU], incOne);
    }
    if(beta != d_one){
        cblas_dscal(_M * _Crk, (beta) , _CU, incOne);
    }
    //hc_printmat(_CU, _M, _M, ld_CU);
    int CU_nrows = _M;
    /*printf("__qra: CU_ncols:%d new_CU_ncols:%d\n", CU_ncols, (CU_ncols+_Ark));*/
    CU_ncols += _Ark;
    *pnew_CU_ncols = CU_ncols; //CHANGED VALUE, RETURN VALUE
    if(gemm_print_index){
        printf(" QRA\t|%d\t|CU_nrows:%d CU_ncols:%d ld_CU:%d  QRA:%p\n", __LINE__, CU_nrows, CU_ncols, ld_CU, _CU);
    }
    info = LAPACKE_dgeqrf(
            LAPACK_COL_MAJOR, CU_nrows, CU_ncols, _CU, ld_CU, qrtauA);
    //printf("%d %d %d _M:%d\n", CU_nrows, CU_ncols, ld_CU, _M);
    //hc_printmat(_CU, _M, _M, ld_CU);
    /*hc_printmat(qrtauA, 1, _M, 1); */
    if(info != 0){
        fprintf(stderr,
                "%s %d ERROR in LAPACKE_dgeqrf(1:CU_nrows:%d 2:CU_ncols:%d 3:_CU:%p 4:ld_CU:%d 5:qrtauA:%p) info=%d maxrank:%d\n",
                __FILE__, __LINE__, CU_nrows, CU_ncols, _CU, ld_CU, qrtauA, info, maxrank);
        exit(-1);
    }
    //hc_printmat(_AU, _M, _M, ld_AU);
}
/*
 * CV, AV and BV are stored as transposed.
 * Performs QR([CV|(AV^T*BU*BV^T)])
 * Step 1: F=AV^T*BU gemm
 * Step 2: G=P*BV^T  gemm
 * Step 3: CV=CV|G   done in the gemm at Step 2
 * Step 4: QR(CV)    potrf
 */
void __qrb(
        int _M,
        int maxrank,
        double* _CV, int ld_CV, int _Crk,
        int* pnew_CV_ncols,
        double* _AV, int ld_AV, int _Ark,
        double* _BU, int ld_BU,
        double* _BV, int ld_BV, int _Brk,
        double* qrtauB, double* AcolBcolT){
    int info;
    /// B2 = ( (iA{2} * iB{2}') * iB{1}');
    assert(AcolBcolT != NULL);
    int ld_AcolBcolT = maxrank; //ASSUMPTION

    double alpha = 1.0;
    double beta = 0.0;
    int AV_ncols = _Ark; //ASSUMPTION
    int BV_nrows = _M;   //ASSUMPTION
    int BV_ncols = _Brk; //ASSUMPTION
    // AV*BV^T
    if(gemm_print_index){
        printf(" QRB\t|%d\t| (AV*BV^T)       M,N,K:%d,%d,%d  LDA,LDB,LDC:%d,%d,%d alpha:%g beta:%g\n", __LINE__,AV_ncols, BV_ncols, BV_nrows,ld_AV,ld_BV,ld_AcolBcolT, alpha, beta);
    }
    /* Step 1: F=AV^T*BV gemm */
    cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
            AV_ncols, BV_ncols, BV_nrows, (alpha),
            _AV, ld_AV, _BV, ld_BV, (beta),
            AcolBcolT, ld_AcolBcolT);
    int AcolBcolT_nrows = AV_ncols;
    int AcolBcolT_ncols = BV_ncols;

    int BU_nrows = _M;     //ASSUMPTION
    int BU_ncols = _Brk;  //ASSUMPTION
    int CV_ncols = _Crk;
    if ((AcolBcolT_nrows + _Crk) > 2*maxrank){
        fprintf(stderr, "%s %s %d: Sum of two ranks (%d) is too big! \
                AcolBcolT:%d _Crk:%d maxrank:%d (x2: %d)\n",
                __FILE__, __func__, __LINE__, (AcolBcolT_nrows + _Crk), AcolBcolT, _Crk, maxrank, 2*maxrank);
        exit(-1);
    }
    // (AV*BV^T) * BU^T
    if(gemm_print_index){
        printf(" QRB\t|%d\t| (AV*BV^T) * BU^T  M,N,K:%d,%d,%d  LDA,LDB,LDC:%d,%d,%d alpha:%g beta:%g   CV_ncols:%d AcolBcolT_nrows:%d ldcB:%d\n",
                __LINE__, BU_nrows,AcolBcolT_nrows, BU_ncols,ld_BU,ld_AcolBcolT,ld_CV, alpha, beta, CV_ncols, AcolBcolT_nrows, ld_CV);
    }
    /* Step 2: G=(F*BU^T)^T  gemm */
    /* Step 3: CV=CV | G   done in the gemm at Step 2*/
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans,
            BU_nrows, AcolBcolT_nrows, BU_ncols,
            (alpha),
            _BU, ld_BU,
            AcolBcolT, ld_AcolBcolT,
            (beta),
            &_CV[CV_ncols*ld_CV], ld_CV);
    CV_ncols        += AcolBcolT_nrows;
    *pnew_CV_ncols   = CV_ncols;       //CHANGED VALUE, RETURN VALUE
    int CV_nrows = _M;             //ASSUMPTION

    if(gemm_print_index){
        printf(" QRB\t|%d\t|CV_nrows:%d CV_ncols:%d ld_CV:%d QRB:%p\n",
                __LINE__, CV_nrows, CV_ncols, ld_CV, _CV);
    }
    /* Step 4: QR(CV)    potrf */
    info = LAPACKE_dgeqrf(LAPACK_COL_MAJOR, CV_nrows, CV_ncols,
            _CV, ld_CV, qrtauB);
    if(info != 0){
        fprintf(stderr,
                "%s %d ERROR in LAPACKE_dgeqrf(1:CV_nrows:%d 2:CV_ncols:%d 3:_CV:%p 4:ld_CV:%d 5:qrtauB:%p) info=%d maxrank:%d\n",
                __FILE__, __LINE__, CV_nrows, CV_ncols, _CV, ld_CV, qrtauB,info,  maxrank);
        exit(-1);
    }
}
void __qrb_b_dense(
        int _M,
        int maxrank,
        double* _CV, int ld_CV, int _Crk,
        int* pnew_CV_ncols,
        double* _AV, int ld_AV, int _Ark,
        double* _B,  int ld_B,
        double* qrtauB){
    int info;
    int CV_ncols = _Crk;

    double alpha = 1.0;
    double beta = 0.0;
    int AV_ncols = _Ark; //ASSUMPTION
    // (AV*BV^T) * BU^T
    if(gemm_print_index){
        printf(" QRB\t|%d\t| AV*B^T  M,N,K:%d,%d,%d  LDA,LDB,LDC:%d,%d,%d alpha:%g beta:%g   CV_ncols:%d ldcB:%d\n",
                __LINE__, _Ark, _M, _M, ld_AV, ld_B, ld_CV, alpha, beta, CV_ncols, ld_CV);
    }
    /* Step 2: G=P*BV^T  gemm */
    /* Step 3: CV=CV | G   done in the gemm at Step 2*/

    /* trans(AV) * trans(B) */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
            _M, 
            _Ark, 
            _M,
            (alpha),
            _B, ld_B,
            _AV, ld_AV,
            (beta),
            &_CV[CV_ncols*ld_CV], ld_CV);
    CV_ncols        += _Ark;
    *pnew_CV_ncols   = CV_ncols;       //CHANGED VALUE, RETURN VALUE
    int CV_nrows = _M;             //ASSUMPTION

    if(gemm_print_index){
        printf(" QRB\t|%d\t|CV_nrows:%d CV_ncols:%d ld_CV:%d QRB:%p\n",
                __LINE__, CV_nrows, CV_ncols, ld_CV, _CV);
    }
    /* Step 4: QR(CV)    potrf */
    info = LAPACKE_dgeqrf(LAPACK_COL_MAJOR, CV_nrows, CV_ncols,
            _CV, ld_CV, qrtauB);
    if(info != 0){
        fprintf(stderr,
                "%s %d ERROR in LAPACKE_dgeqrf(1:CV_nrows:%d 2:CV_ncols:%d 3:_CV:%p 4:ld_CV:%d 5:qrtauB:%p) info=%d maxrank:%d\n",
                __FILE__, __LINE__, CV_nrows, CV_ncols, _CV, ld_CV, qrtauB,info,  maxrank);
        exit(-1);
    }
}

//void __svd(ka_matrix* _CU, ka_matrix* _CV, int rank, double acc, ka_matrix* _U, ka_matrix* _V)
// Uses CV for TRMM
void __svd(
        int _M,
        int maxrank,/*this is compmaxrank*/
        double* _CU, int ld_CU,
        double* _CV, int ld_CV, int _Crk,
        double* _U, int ld_U,
        double* _V, int ld_V, int* pnew_UVrk,
        int rank,
        double acc,
        double* _rA, double* _rB, double* _T, double* sigma, double* svdsuperb
        ) {
    int info;
    int nb        = _M;   //ASSUMPTION
    int CU_nrows  = _M;   //ASSUMPTION
    int CU_ncols  = _Crk; //ASSUMPTION
    /// rA     START
    int rA_nrows  = CU_nrows < CU_ncols ? CU_nrows : CU_ncols;
    int rA_ncols  = CU_ncols;    //ASSUMPTION
    int maxncolsR = 2*_Crk; //nb;          //ASSUMPTION
    int ld_rA     = rA_nrows;

    if(rA_nrows != rA_ncols){
        printf("TRMM cannot be used because R coming from QR factorization of A is not square nrows: %d ncols:%d \n", rA_nrows, rA_ncols);
        exit(1);
    }
    assert(_rA != NULL);

    if(gemm_print_index){
        printf(" SVD\t|%d\t| copy rA rA_nrows:%d rA_ncols:%d ld_CU:%d ld_rA:%d CU:%p rA:%p\n",
              __LINE__, rA_nrows, rA_ncols, ld_CU, ld_rA, _CU, _rA);
    }
    double zero = 0.0;
    char chlow = 'L';
    LAPACK_dlaset(&chlow,
            &rA_nrows, &rA_ncols, &zero, &zero, _rA, &ld_rA);
    char chup = 'U';
    LAPACK_dlacpy(&chup,
            &rA_nrows, &rA_ncols,
            _CU, &ld_CU,
            _rA, &ld_rA);
    if(gemm_print_mat){
        printf("%d\t|_CU and _rA\n", __LINE__);
        hc_printmat(_CU,  _M, _Crk, ld_CU);
        hc_printmat(_rA,  rA_nrows, rA_ncols, ld_rA);
    }
    /// rA     END
    /// rB    START
    int CV_nrows = _M;         //ASSUMPTION
    int CV_ncols = _Crk;       //ASSUMPTION
    int rB_nrows = CV_nrows < CV_ncols ? CV_nrows : CV_ncols;
    int rB_ncols = CV_ncols;
    int ld_rB    = rB_nrows;     //ASSUMPTION
    assert(rA_ncols == rB_ncols);
    // trmm does not access lower part but gemm reads.
    // However lower part is never written.
    if(gemm_print_index){
        printf(" SVD\t|%d\t| copy rB rB_nrows:%d rB_ncols:%d ld_rB:%d ld_CV:%d rB:%p CV:%p\n",
              __LINE__, rB_nrows, rB_ncols, ld_rB, ld_CV, _rB, _CV);
    }
    if(use_trmm == 0){
        assert(_rB != NULL);
        LAPACK_dlaset(&chlow,
                &rB_nrows, &rB_ncols, &zero, &zero, _rB, &ld_rB);
        LAPACK_dlacpy(&chup,
                &rB_nrows, &rB_ncols, _CV, &ld_CV, _rB, &ld_rB);
    } else {
        _rB = _CV;
        ld_rB = ld_CV;
    }
    if(gemm_print_mat){
        printf("%d\t|_CV and _rB\n", __LINE__);
        hc_printmat(_CV,  _M, _Crk, ld_CV);
        hc_printmat(_rB,  rB_nrows, rB_ncols, ld_rB);
    }
    /// rB    END
    /// SVD   START
    int T_nrows = rA_ncols;
    int T_ncols = rA_ncols;
    int ld_T    = T_nrows;
    double alpha = 1.0;
    double beta  = 0.0;
    if(use_trmm == 1){
        cblas_dtrmm(CblasColMajor, CblasRight, CblasUpper, CblasTrans, CblasNonUnit, rA_nrows, rA_ncols,  alpha, _rB, ld_rB, _rA, ld_rA); //FIXME Correctness of rA_ncols is not checked
        _T = _rA;
    } else {
        assert(_T != NULL);
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans,
                rA_nrows, rB_nrows, rB_ncols, (alpha), _rA, ld_rA,
                _rB, ld_rB, (beta), _T, ld_T);
    }
    if(gemm_print_index){
        printf(" SVD\t|%d\t| T=rA*rB^T  rA_nrows:%d rB_nrows:%d rB_ncols:%d ld_rA:%d ld_rB:%d ld_T:%d alpha:%g beta:%g\n",
              __LINE__, rA_nrows, rB_nrows, rB_ncols, ld_rA, ld_rB, ld_T, alpha, beta);
    }
    if(gemm_print_mat){
        hc_printmat(_T,  T_nrows, T_ncols, ld_T);
    }
    // Singular values are ALWAYS floating point numbers.
    assert(sigma != NULL);
    int size_sigma = T_nrows;
    assert(svdsuperb != NULL);
    if(gemm_print_index){
        printf(" SVD\t|%d\t| svd(T)    (3.m)T_nrows:%d (4.n)T_ncols:%d ld_T:%d ld_U:%d (11.ldvt)ld_V:%d _T:%p (zero based parameter indices)\n",
              __LINE__, T_nrows, T_ncols,  ld_T, ld_U, ld_V, _T);
    }
    info = LAPACKE_dgesvd(LAPACK_COL_MAJOR, 'A', 'A',
            T_nrows, T_ncols, _T, ld_T, sigma,
            _U, ld_U, _V, ld_V,
            svdsuperb);
    if(info != 0){
        fprintf(stderr,
                "%s %d ERROR in LAPACKE_dgesvd() info=%d "
                "1:T_nrows=%d, 2:T_ncols=%d, 3:_T=%p, 4:ld_T=:%d, 5:sigma=%p, "
                "6:_U=%p, 7:ld_U=%d, 8:_V=%p, 9:ld_V:%d, "
                "10:svdsuperb:%p "
                "\n",
                __FILE__, __LINE__, info,
                T_nrows, T_ncols, _T, ld_T, sigma,
                _U, ld_U, _V, ld_V,
                svdsuperb);
        exit(-1);
    }
    int U_nrows, U_ncols, V_nrows, V_ncols;
    U_nrows = U_ncols = V_nrows = V_ncols = T_nrows;

    if(gemm_print_mat){
        hc_printmat(_U,  U_nrows, U_ncols, ld_U);
        hc_printmat(sigma,  1, size_sigma, 1);
        hc_printmat(_V,  V_nrows, V_ncols, ld_V);
        printf("%d %e\n", rank, acc);
    }
    int finalrank = -1;
    //double relacc = (acc*sigma[0]*acc);
    double relacc = (acc);
    //double relacc = (acc*acc);
    if(rank != 0) { /// truncate according to rank
        finalrank = rank;
        if(rank > size_sigma)
            finalrank = size_sigma;
    }
    else{ /// truncate according to acc
        int newrank = size_sigma;
        int i;
        for(i=2;i<size_sigma;i++){
            if(sigma[i] < relacc)
            {
                newrank=i;
                break;
            }
        }
        finalrank = newrank;
    }
    if(finalrank > maxrank){
        fprintf(stderr, "%s %s %d: Rank after truncation is too big! finalrank:%d maxrank:%d. Increase computation maxrank.\n", __FILE__, __func__, __LINE__, finalrank, maxrank);
        exit(-1);
    }

    if(gemm_print_index){int i;
        printf("rank:%d acc:%.2e   relac:%.2e size_sigma:%d final_rank:%d:  ",
                rank,   acc, relacc, size_sigma, finalrank);
        for(i=0;i<size_sigma;i++){
            printf("%d:%.2e ", i,sigma[i]);
        }
        printf("\n");
    }
    if(gemm_print_index){
        printf("size_sigma:%d finalrank:%d %.2e\n", size_sigma, finalrank, acc);
    }
    U_ncols = finalrank;
    V_nrows = finalrank;
    /// SVD    END
    int rank_V = finalrank;
    int k;
    for(k = 0; k < rank_V; k++){
        double diagval = sigma[k];
        cblas_dscal(V_ncols, (diagval), &_V[k], ld_V);
    }
    if(gemm_print_index){
    printf(" SVD\t|%d\t| S*V     V_ncols:%d ld_V:%d _V:%p\n",
              __LINE__, V_ncols,  ld_V, _V);
    }
    if(gemm_print_mat){
        hc_printmat(_V,  V_nrows, V_ncols, ld_V);
    }
    *pnew_UVrk = finalrank;
}
//void __newu(ka_matrix* _CU, double* qrtauA, ka_matrix* _U)
void __newu(
        int _M,
        int ncols_qA,
        double* _CU, int ld_CU, int _Crk,
        double* _U,  int ld_U,  int _Urk,
        double* qrtauA
        ) {
    int info = 0;
    int CU_nrows = _M;
    int CU_ncols = ncols_qA;
    int U_nrows = ncols_qA; //ASSUMPTION
    int U_ncols = _Urk; //ASSUMPTION
    int nrows = CU_nrows - U_nrows;
    double zero  = 0.0;
    if(gemm_print_index){
        printf(" NEWU\t|%d\t|    zero     nrows:%d U_ncols:%d ld_U:%d _Crk:%d CU_ncols:%d _Urk:%d   CU_nrows:%d U_nrows:%d diff:%d\n",
              __LINE__, nrows, U_ncols,  ld_U, _Crk, CU_ncols,  _Urk, CU_nrows, U_nrows, nrows);
    }
    //info = LAPACKE_zlaset(LAPACK_COL_MAJOR, 'A', nrows, U_ncols,
    //        zero, zero, &_U[U_nrows], ld_U);

    char uplo = 'A';
    LAPACK_dlaset( &uplo, &nrows, &U_ncols, &zero, &zero, &_U[U_nrows], &ld_U );
    if(info != 0){
        fprintf(stderr,
                "%s %d ERROR in LAPACKE_dlaset() info=%d\n",
                __FILE__, __LINE__, info);
    }
    if(gemm_print_mat){
        hc_printmat(_U,  _M, _M, ld_U);
    }

    //zunmqr
    info = LAPACKE_dormqr(LAPACK_COL_MAJOR, 'L', 'N',
            CU_nrows, U_ncols, ncols_qA, _CU, ld_CU, qrtauA, _U, ld_U);
    if(gemm_print_index){
        printf(" NEWU\t|%d\t|    ormqr     CU_nrows (new U_nrows):%d U_ncols:%d ncols_qA:%d ld_CU:%d ld_U:%d\n",
              __LINE__, CU_nrows, U_ncols, ncols_qA, ld_CU, ld_U);
    }
    if(gemm_print_mat){
        hc_printmat(_U,  _M, _M, ld_U);
    }
    if(info != 0){
        fprintf(stderr,
                "%s %d ERROR in LAPACKE_dormqr() info=%d\n",
                __FILE__, __LINE__, info);
        printf(" NEWU\t|%d\t|    ormqr     CU_nrows (new U_nrows):%d U_ncols:%d ncols_qA:%d ld_CU:%d ld_U:%d U_nrows:%d\n",
              __LINE__, CU_nrows, U_ncols, ncols_qA, ld_CU, ld_U, U_nrows);
        if(0){
            int i, j, ssend, ldarr;
            double *arr;
            arr = _CU; ldarr = ld_CU; ssend = ncols_qA;
            for(i=0;i<50;i++){
                for(j=0;j<3;j++){
                    printf("%.3e ", arr[j*ldarr+i]);
                }
                printf("...");
                for(j=ssend-4;j<ssend;j++){
                    printf("%.3e ", arr[j*ldarr+i]);
                }
                printf("\n");
            }
            printf("\n");
            arr = _U; ldarr = ld_U; ssend = U_ncols;
            for(i=0;i<50;i++){
                for(j=0;j<3;j++){
                    printf("%.3e ", arr[j*ldarr+i]);
                }
                printf("...");
                for(j=ssend-4;j<ssend;j++){
                    printf("%.3e ", arr[j*ldarr+i]);
                }
                printf("\n");
            }
            for(j=0;j<U_ncols;j++){
                for(i=0; i < CU_nrows; i++){
                    double val = _U[j*ld_U+i];
                    if(val != val){
                        printf("%d,%d is nan (%g) CU_nrows:%d U_ncols:%d ld_U:%d\n", i, j, val,  CU_nrows, U_ncols, ld_U);
                    }
                }
            }
        }
        exit(-1);
    }
    U_nrows = CU_nrows;
    LAPACKE_dlacpy(LAPACK_COL_MAJOR, 'A', U_nrows, U_ncols,
            _U, ld_U, _CU, ld_CU);
    if(gemm_print_index){
        printf(" NEWU\t|%d\t|    copy     U_nrows:%d U_ncols:%d ld_CU:%d ld_U:%d\n",
                __LINE__, U_nrows, U_ncols, ld_CU, ld_U);
        if(0){
            int i, j;
            for(i=0;i<100;i++){
                for(j=0;j<3;j++){
                    printf("%.3e ", _CU[j*ld_CU+i]);
                }
                printf("...");
                for(j=U_ncols-4;j<U_ncols;j++){
                    printf("%.3e ", _CU[j*ld_CU+i]);
                }
                printf("\n");
            }
            getc(stdin);
        }
    }
}
//void __newv(ka_matrix* _CV, double* qrtauB, ka_matrix* _V)
void __newv(
        int _M,
        int ncols_qB,
        double* _CV, int ld_CV, int _Crk,
        double* _V,  int ld_V,  int _Vrk,
        double* qrtauB
        ) {
    int info;
    int CV_nrows = _M;
    int CV_ncols = ncols_qB;
    int V_nrows = _Vrk; //ASSUMPTION
    int V_ncols = ncols_qB; //ASSUMPTION
    int ncols = CV_nrows - V_ncols;
    double zero  = 0.0;
    if(gemm_print_index){
        printf(" NEWV\t|%d\t|    zero     V_nrows:%d ncols:%d ld_V:%d _Crk:%d CV_ncols:%d _Vrk:%d CV_nrows:%d V_ncols:%d diff:%d\n",
              __LINE__, V_nrows,  ncols, ld_V, _Crk, CV_ncols,  _Vrk, CV_nrows, V_ncols, ncols);
    }
    //LAPACKE_zlaset(LAPACK_COL_MAJOR, 'A', V_nrows, ncols,
    //        zero, zero, &(_V[V_ncols*ld_V]), ld_V);
    char uplo = 'A';
    size_t iv = V_ncols*ld_V;
    LAPACK_dlaset( &uplo, &V_nrows, &ncols, &zero, &zero, &(_V[iv]), &ld_V );
    if(gemm_print_mat){
        hc_printmat(_V,  _M, _M, ld_V);
    }
    if(gemm_print_index){
        printf(" NEWV\t|%d\t|    ormqr     V_nrows:%d CV_nrows:%d ncols_qB:%d ld_CV:%d ld_V:%d\n",
              __LINE__, V_nrows, CV_nrows, ncols_qB, ld_CV, ld_V);
    }
    //zunmqr
    info = LAPACKE_dormqr(LAPACK_COL_MAJOR, 'R', 'T',
            V_nrows, CV_nrows, ncols_qB, _CV, ld_CV, qrtauB, _V, ld_V);
    if(gemm_print_mat){
        hc_printmat(_V,  _M, _M, ld_V);
    }
    if(info != 0){
        fprintf(stderr,
                "%s %d ERROR in LAPACKE_dormqr() info=%d\n",
                __FILE__, __LINE__, info);
        exit(-1);
    }
    V_ncols  = CV_nrows;
    CV_ncols = V_nrows;
    if(gemm_print_index){
        printf(" NEWV\t|%d\t|    trans    V_nrows:%d V_ncols:%d ld_V:%d ld_CV%d\n",
              __LINE__, V_nrows,  V_ncols, ld_V, ld_CV);
    }
    LAPACKE_dge_trans(LAPACK_COL_MAJOR, V_nrows, V_ncols,
            _V, ld_V, _CV, ld_CV);
    if(gemm_print_index){
        printf(" NEWV\t|%d\t|    copy     V_nrows:%d V_ncols:%d ld_CV:%d ld_V:%d\n",
              __LINE__, V_nrows, V_ncols, ld_CV, ld_V);
          if(0){
              int i, j;
              for(i=0;i<100;i++){
                  for(j=0;j<3;j++){
                      printf("%.3e ", _CV[j*ld_CV+i]);
                  }
                  printf("...");
                  for(j=V_ncols-4;j<V_ncols;j++){
                      printf("%.3e ", _CV[j*ld_CV+i]);
                  }
                  printf("\n");
              }
            getc(stdin);
          }
    }
    if(gemm_print_mat){
        hc_printmat(_CV,  _M, _M, ld_CV);
    }
}
/**
 * Performs C = alpha A B + beta C
 * Each of A, B and C is stored as U*V^T.
 * alpha and beta are scalars.
 *
 */
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
        int storage_maxrank,
        int maxrk, /*this is compmaxrank*/
        double acc,
        double* work
)
{
    if(gemm_print_index){
		printf("%d:%s work:%p ", __LINE__, __func__, work);
		printf("M:%d N:%d  LDA:%d LDB:%d LDC:%d rk:%d maxrk:%d acc:%e a:%e b:%e\n",
				M, N, LDA, LDB, LDC, rk, maxrk, acc, alpha, beta);
	}

    int new_Crk = 0;
    /*
     * NOTES:
     * assumptions on matrix dimensions are marked as //ASSUMPTION
     */
    /*printf("%d %d|%g->%d. %g %g, %g %g\n",  */
            /*__LINE__, __COUNTER__,*/
            /*Crk[0], new_Crk, CU[0], CU[1], CV[0], CV[1]);*/
    //hcore_dgemm(Aij, ik, Ajk, -1, rank, acc);
    int _Ark = (int)(Ark[0]);
    int _Brk = (int)(Brk[0]);
    int _Crk = (int)(Crk[0]);
    if(_Ark == 0 || _Brk == 0 || _Crk == 0){
        fprintf(stderr, "%s %d: _Ark=%d _Brk=%d _Crk=%d. These rank values should not be zero.\n", __FILE__, __LINE__, _Ark, _Brk, _Crk);
        exit(-1);
        //return MORSE_ERR_ILLEGAL_VALUE; 
    }

    int _M = M; int _N = N;
    double* _CU = CU; int ld_CU = LDC;
    double* _CV = CV; int ld_CV = LDC; int* pnew_Crk = &new_Crk;
    double* _AU = AU; int ld_AU = LDA;
    double* _AV = AV; int ld_AV = LDA;
    double* _BU = BU; int ld_BU = LDB;
    double* _BV = BV; int ld_BV = LDB;
    int rank = rk;


    //printf("%s %s %d maxrank=%d\n", __FILE__, __func__, __LINE__, maxrk);
    char chall = 'A';
    double* CUclone = NULL;
    size_t CUclone_nelm =  _M * 2 * maxrk;
    int ld_CUclone = _M;
    int use_CUV_clone = 1;
    if(use_CUV_clone == 1) {
        if(use_scratch){
            CUclone = work;
            work += CUclone_nelm;
        } else {
            CUclone = malloc(CUclone_nelm * sizeof(double));
        }
        LAPACK_dlacpy(&chall,
                &_M, &_Crk,
                _CU, &ld_CU,
                CUclone, &ld_CUclone);
    }

    double* CVclone = NULL;
    size_t CVclone_nelm = _M * 2 * maxrk;
    int ld_CVclone = _M;
    if(use_CUV_clone == 1) {
        if(use_scratch){
            CVclone = work;
            work += CVclone_nelm;
        } else {
            CVclone = malloc(CVclone_nelm * sizeof(double));
        }
        LAPACK_dlacpy(&chall,
                &_M, &_Crk,
                _CV, &ld_CV,
                CVclone, &ld_CVclone);
    }
    double* _CU_save = _CU;
    double* _CV_save = _CV;

    if(use_CUV_clone == 1) {
        _CU = CUclone;
        _CV = CVclone;
    }
    int nb    = _M;//ASSUMPTION
    double* qrtauA = NULL;
    size_t qrtauA_nelm = nb;
    if(use_scratch){
        qrtauA = work;
    } else {
        qrtauA = malloc(qrtauA_nelm * sizeof(double));
    }
    assert(qrtauA != NULL);
    double* qrtauB = NULL;
    size_t qrtauB_nelm = nb;
    if(use_scratch){
        qrtauB = work + qrtauA_nelm;
    } else {
        qrtauB = malloc(qrtauB_nelm * sizeof(double));
    }
    assert(qrtauB != NULL);



    int CU_ncols = 0;
    //Warning:In this function, there are assumptions
    //on leading dimensions and number of rows/cols
    __qra(_M, maxrk, _CU, ld_CU, _Crk, &CU_ncols, _AU, ld_AU, _Ark, alpha, beta, qrtauA);
    //output: _CU contains [_CU _AU]. use _Crk+_Ark as number of cols
    assert(CU_ncols == (_Crk + _Ark));

    int CV_ncols = 0;
    //Warning:In this function, there are assumptions
    //on leading dimensions and number of rows/cols
    double* qrb_aubut = NULL;
    size_t qrb_aubut_nelm = maxrk * maxrk;
    if(use_scratch){
       qrb_aubut = work + qrtauA_nelm + qrtauB_nelm;
    } else {
       qrb_aubut = malloc(qrb_aubut_nelm * sizeof(double));
    }
    __qrb(_M, maxrk, _CV, ld_CV, _Crk, &CV_ncols,
            _AV, ld_AV, _Ark, _BU, ld_BU, _BV, ld_BV, _Brk, qrtauB, qrb_aubut);
    if(CU_ncols == 0 || CV_ncols == 0){
        fprintf(stderr, "%s %d: CU_ncols=%d CV_ncols=%d. These values should not be zero.\n", __FILE__, __LINE__, CU_ncols, CV_ncols);
        exit(-1);
        //return MORSE_ERR_ILLEGAL_VALUE; 
    }
    if(use_scratch == 0) {
        free(qrb_aubut);
    }
    assert(CU_ncols == CV_ncols);


    int ld_newU = nb; //ASSUMPTION
    int ld_newV = 2*maxrk; //ASSUMPTION
    int new_UVrk;
    double* newU = NULL;
    size_t newU_nelm = ld_newU * 2*maxrk;
    if(use_scratch){
        newU = work + qrtauA_nelm + qrtauB_nelm + qrb_aubut_nelm;
    } else {
        newU = malloc(newU_nelm * sizeof(double));
    }
    double* newV = NULL;
    size_t newV_nelm = nb * ld_newV;

    if(use_scratch){
        newV = work + qrtauA_nelm + qrtauB_nelm + qrb_aubut_nelm + newU_nelm;
    } else {
        newV = malloc(newV_nelm *  sizeof(double));
    }
    assert(newU != NULL);
    assert(newV != NULL);
    double *svd_rA = NULL;
    int svd_rA_nrows  = _M < CU_ncols ? _M : CU_ncols;
    int svd_rA_ncols  = CU_ncols;    //ASSUMPTION
    size_t svd_rA_nelm;
    if(use_trmm == 1){ // allocate more because rA will be output of trmm(rA*rB^T)
        svd_rA_nelm = svd_rA_nrows * svd_rA_nrows;
    } else {
        svd_rA_nelm = svd_rA_nrows * svd_rA_ncols;
    }
    if(use_scratch == 1) {
        svd_rA = work + qrtauA_nelm + qrtauB_nelm + qrb_aubut_nelm + newU_nelm + newV_nelm;
    } else {
        svd_rA = malloc(svd_rA_nelm * sizeof(double));
    }
    double *svd_rB = NULL;
    int svd_rB_nrows = _M < CV_ncols ? _M : CV_ncols;
    int svd_rB_ncols = CV_ncols;
    size_t svd_rB_nelm;
    if(use_trmm == 0) {
        svd_rB_nelm = svd_rB_nrows * svd_rB_ncols;
        if(use_scratch == 1) {
            svd_rB = work + qrtauA_nelm + qrtauB_nelm + qrb_aubut_nelm + newU_nelm + newV_nelm + svd_rA_nelm;
        } else {
            svd_rB = malloc(svd_rB_nelm * sizeof(double));
        }
    } else {
        svd_rB_nelm = 0; //This variable will be used for calculating offset for work array
    }
    double *svd_T = NULL;
    int svd_T_nrows = svd_rA_ncols;
    int svd_T_ncols = svd_rA_ncols;
    size_t svd_T_nelm;
    if(use_trmm == 0) {
        svd_T_nelm = svd_T_nrows * svd_T_ncols;
        if(use_scratch == 1) {
            svd_T = work + qrtauA_nelm + qrtauB_nelm + qrb_aubut_nelm + newU_nelm + newV_nelm + svd_rA_nelm + svd_rB_nelm;
        } else {
            svd_T = malloc(svd_T_nelm * sizeof(double));
        }
    } else {
        svd_T_nelm = 0; //This variable will be used for calculating offset for work array
    }

    double *svd_sigma  = NULL;
    double *svd_superb = NULL;
    size_t svd_sigma_nelm  = svd_T_nrows;
    size_t svd_superb_nelm = svd_T_nrows;
    if(use_scratch == 1){
        svd_sigma  = work + qrtauA_nelm + qrtauB_nelm + qrb_aubut_nelm + newU_nelm + newV_nelm + svd_rA_nelm + svd_rB_nelm + svd_T_nelm;
        svd_superb = work + qrtauA_nelm + qrtauB_nelm + qrb_aubut_nelm + newU_nelm + newV_nelm + svd_rA_nelm + svd_rB_nelm + svd_T_nelm + svd_sigma_nelm;
    } else {
        svd_sigma  = malloc(svd_sigma_nelm  * sizeof(double));
        svd_superb = malloc(svd_superb_nelm * sizeof(double));
    }
    if(ld_newV < CU_ncols){ /* ld_newV=2*maxrank CU_ncols=Crk+Ark */
        fprintf(stderr, "%s %d: Increase computation maxrank. %d is not enough ld_newV:%d CU_ncols:%d\n", __FILE__, __LINE__, maxrk, ld_newV, CU_ncols);
        exit(-1);
        //return MORSE_ERR_ILLEGAL_VALUE; 
    }
    __svd(
            _M,
            maxrk,
            _CU,  ld_CU,
            _CV,  ld_CV, CU_ncols,
            newU,   ld_newU,
            newV,   ld_newV,  &new_UVrk ,
            rank, acc,
            svd_rA, svd_rB, svd_T, svd_sigma, svd_superb
         );
    if(use_scratch == 0) {
        free(svd_rA);
    }
    if(use_scratch == 0 && use_trmm == 0){
        free(svd_rB);
    //    free(svd_T);
    }
    int ncols_qA = CU_ncols;
    //__newu(_CU, qrtauA, newU);
    if(new_UVrk > storage_maxrank) {
        fprintf(stderr, "%s %d: Increase storage maxrank %d. Observed rank for U and V is %d\n", __FILE__, __LINE__, storage_maxrank, new_UVrk);
        exit(-1);
    }
    __newu(_M,
            ncols_qA,
            _CU, ld_CU, _Crk,
            newU, ld_newU, new_UVrk,
            qrtauA
            );
    //Warning: number of columns of _CU is new_UVrk now!

    int ncols_qB = CV_ncols;
    //__newv(_CV, qrtauB, newV);
    __newv(_M,
            ncols_qB,
            _CV, ld_CV, _Crk,
            newV, ld_newV, new_UVrk,
            qrtauB
            );
    *pnew_Crk = new_UVrk;
    //Warning: number of columns of _CV is new_UVrk now!

    //printf("%d->%d\n", _Crk, *pnew_Crk);

    if(use_CUV_clone == 1) {
        LAPACK_dlacpy(&chall,
                &_M, &new_UVrk,
                CUclone, &ld_CUclone,
                _CU_save, &ld_CU
                );

        LAPACK_dlacpy(&chall,
                &_M, &new_UVrk,
                CVclone, &ld_CVclone,
                _CV_save, &ld_CV
                );
        if(use_scratch == 0) {
            free(CUclone);
            free(CVclone);
        }
    }


    if(use_scratch == 0) {
        free(qrtauA);
        free(qrtauB);
        free(newU);
        free(newV);
        free(svd_sigma);
        free(svd_superb);
    }


    int old_Crk = Crk[0];
    if(gemm_print_index){
        printf("Ark:%d Brk:%d Crk[0]:%d %g RANK CHANGE: %d->%d\n",
                _Ark, _Brk, Crk[0], Crk[0], old_Crk, *pnew_Crk);
    }
    Crk[0] = new_Crk;
    if(gemm_print_index){
        int casted_Crk = (int)(Crk[0]);
        printf("casted_Crk:%d Ark:%d Brk:%d Crk[0]:%d %g RANK CHANGE: %d->%d\n",
                casted_Crk, _Ark, _Brk, Crk[0], Crk[0], old_Crk, new_Crk);
    }
}
/**
 * Performs C = alpha A B + beta C
 * Each of A, B and C is stored as U*V^T.
 * alpha and beta are scalars.
 *
 */
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
)
{
    if(gemm_print_index){
		printf("%d:%s work:%p ", __LINE__, __func__, work);
		printf("M:%d N:%d  LDA:%d LDB:%d LDC:%d rk:%d maxrk:%d acc:%e a:%e b:%e\n",
				M, N, LDA, LDB, LDC, rk, maxrk, acc, alpha, beta);
	}

    int new_Crk = 0;
    /*
     * NOTES:
     * assumptions on matrix dimensions are marked as //ASSUMPTION
     */
    /*printf("%d %d|%g->%d. %g %g, %g %g\n",  */
            /*__LINE__, __COUNTER__,*/
            /*Crk[0], new_Crk, CU[0], CU[1], CV[0], CV[1]);*/
    //hcore_dgemm(Aij, ik, Ajk, -1, rank, acc);
    int _Ark = (int)(Ark[0]);
    int _Crk = (int)(Crk[0]);
    if(_Ark == 0 || _Crk == 0){
        fprintf(stderr, "%s %d: _Ark=%d _Crk=%d. These rank values should not be zero.\n", __FILE__, __LINE__, _Ark, _Crk);
        exit(-1);
        //return MORSE_ERR_ILLEGAL_VALUE; 
    }

    int _M = M; int _N = N;
    double* _CU = CU; int ld_CU = LDC;
    double* _CV = CV; int ld_CV = LDC; int* pnew_Crk = &new_Crk;
    double* _AU = AU; int ld_AU = LDA;
    double* _AV = AV; int ld_AV = LDA;
    double* _B  = B;  int ld_B  = LDB;
    int rank = rk;


    //printf("%s %s %d maxrank=%d\n", __FILE__, __func__, __LINE__, maxrk);
    char chall = 'A';
    double* CUclone = NULL;
    size_t CUclone_nelm =  _M * 2 * maxrk;
    int ld_CUclone = _M;
    int use_CUV_clone = 1;
    if(use_CUV_clone == 1) {
        if(use_scratch){
            CUclone = work;
            work += CUclone_nelm;
        } else {
            CUclone = malloc(CUclone_nelm * sizeof(double));
        }
        LAPACK_dlacpy(&chall,
                &_M, &_Crk,
                _CU, &ld_CU,
                CUclone, &ld_CUclone);
    }

    double* CVclone = NULL;
    size_t CVclone_nelm = _M * 2 * maxrk;
    int ld_CVclone = _M;
    if(use_CUV_clone == 1) {
        if(use_scratch){
            CVclone = work;
            work += CVclone_nelm;
        } else {
            CVclone = malloc(CVclone_nelm * sizeof(double));
        }
        LAPACK_dlacpy(&chall,
                &_M, &_Crk,
                _CV, &ld_CV,
                CVclone, &ld_CVclone);
    }
    double* _CU_save = _CU;
    double* _CV_save = _CV;

    if(use_CUV_clone == 1) {
        _CU = CUclone;
        _CV = CVclone;
    }
    int nb    = _M;//ASSUMPTION
    double* qrtauA = NULL;
    size_t qrtauA_nelm = nb;
    if(use_scratch){
        qrtauA = work;
    } else {
        qrtauA = malloc(qrtauA_nelm * sizeof(double));
    }
    assert(qrtauA != NULL);
    double* qrtauB = NULL;
    size_t qrtauB_nelm = nb;
    if(use_scratch){
        qrtauB = work + qrtauA_nelm;
    } else {
        qrtauB = malloc(qrtauB_nelm * sizeof(double));
    }
    assert(qrtauB != NULL);



    int CU_ncols = 0;
    //Warning:In this function, there are assumptions
    //on leading dimensions and number of rows/cols
    __qra(_M, maxrk, _CU, ld_CU, _Crk, &CU_ncols, _AU, ld_AU, _Ark, alpha, beta, qrtauA);
    //output: _CU contains [_CU _AU]. use _Crk+_Ark as number of cols
    assert(CU_ncols == (_Crk + _Ark));

    int CV_ncols = 0;
    //Warning:In this function, there are assumptions
    //on leading dimensions and number of rows/cols
    double* qrb_aubut = NULL;
    size_t qrb_aubut_nelm = maxrk * maxrk;
    if(use_scratch){
       qrb_aubut = work + qrtauA_nelm + qrtauB_nelm;
    } else {
       qrb_aubut = malloc(qrb_aubut_nelm * sizeof(double));
    }
    __qrb_b_dense(_M, maxrk, _CV, ld_CV, _Crk, &CV_ncols,
            _AV, ld_AV, _Ark, _B, ld_B, qrtauB); /* NOT USED, qrb_aubut);*/
    if(CU_ncols == 0 || CV_ncols == 0){
        fprintf(stderr, "%s %d: CU_ncols=%d CV_ncols=%d. These values should not be zero.\n", __FILE__, __LINE__, CU_ncols, CV_ncols);
        exit(-1);
        //return MORSE_ERR_ILLEGAL_VALUE; 
    }
    if(use_scratch == 0) {
        free(qrb_aubut);
    }
    assert(CU_ncols == CV_ncols);


    int ld_newU = nb; //ASSUMPTION
    int ld_newV = 2*maxrk; //ASSUMPTION
    int new_UVrk;
    double* newU = NULL;
    size_t newU_nelm = ld_newU * 2*maxrk;
    if(use_scratch){
        newU = work + qrtauA_nelm + qrtauB_nelm + qrb_aubut_nelm;
    } else {
        newU = malloc(newU_nelm * sizeof(double));
    }
    double* newV = NULL;
    size_t newV_nelm = nb * ld_newV;

    if(use_scratch){
        newV = work + qrtauA_nelm + qrtauB_nelm + qrb_aubut_nelm + newU_nelm;
    } else {
        newV = malloc(newV_nelm *  sizeof(double));
    }
    assert(newU != NULL);
    assert(newV != NULL);
    double *svd_rA = NULL;
    int svd_rA_nrows  = _M < CU_ncols ? _M : CU_ncols;
    int svd_rA_ncols  = CU_ncols;    //ASSUMPTION
    size_t svd_rA_nelm;
    if(use_trmm == 1){ // allocate more because rA will be output of trmm(rA*rB^T)
        svd_rA_nelm = svd_rA_nrows * svd_rA_nrows;
    } else {
        svd_rA_nelm = svd_rA_nrows * svd_rA_ncols;
    }
    if(use_scratch == 1) {
        svd_rA = work + qrtauA_nelm + qrtauB_nelm + qrb_aubut_nelm + newU_nelm + newV_nelm;
    } else {
        svd_rA = malloc(svd_rA_nelm * sizeof(double));
    }
    double *svd_rB = NULL;
    int svd_rB_nrows = _M < CV_ncols ? _M : CV_ncols;
    int svd_rB_ncols = CV_ncols;
    size_t svd_rB_nelm;
    if(use_trmm == 0) {
        svd_rB_nelm = svd_rB_nrows * svd_rB_ncols;
        if(use_scratch == 1) {
            svd_rB = work + qrtauA_nelm + qrtauB_nelm + qrb_aubut_nelm + newU_nelm + newV_nelm + svd_rA_nelm;
        } else {
            svd_rB = malloc(svd_rB_nelm * sizeof(double));
        }
    } else {
        svd_rB_nelm = 0; //This variable will be used for calculating offset for work array
    }
    double *svd_T = NULL;
    int svd_T_nrows = svd_rA_ncols;
    int svd_T_ncols = svd_rA_ncols;
    size_t svd_T_nelm;
    if(use_trmm == 0) {
        svd_T_nelm = svd_T_nrows * svd_T_ncols;
        if(use_scratch == 1) {
            svd_T = work + qrtauA_nelm + qrtauB_nelm + qrb_aubut_nelm + newU_nelm + newV_nelm + svd_rA_nelm + svd_rB_nelm;
        } else {
            svd_T = malloc(svd_T_nelm * sizeof(double));
        }
    } else {
        svd_T_nelm = 0; //This variable will be used for calculating offset for work array
    }

    double *svd_sigma  = NULL;
    double *svd_superb = NULL;
    size_t svd_sigma_nelm  = svd_T_nrows;
    size_t svd_superb_nelm = svd_T_nrows;
    if(use_scratch == 1){
        svd_sigma  = work + qrtauA_nelm + qrtauB_nelm + qrb_aubut_nelm + newU_nelm + newV_nelm + svd_rA_nelm + svd_rB_nelm + svd_T_nelm;
        svd_superb = work + qrtauA_nelm + qrtauB_nelm + qrb_aubut_nelm + newU_nelm + newV_nelm + svd_rA_nelm + svd_rB_nelm + svd_T_nelm + svd_sigma_nelm;
    } else {
        svd_sigma  = malloc(svd_sigma_nelm  * sizeof(double));
        svd_superb = malloc(svd_superb_nelm * sizeof(double));
    }
    if(ld_newV < CU_ncols){ /* ld_newV=2*maxrank CU_ncols=Crk+Ark */
        fprintf(stderr, "%s %d: Increase computation maxrank. %d is not enough ld_newV:%d CU_ncols:%d\n", __FILE__, __LINE__, maxrk, ld_newV, CU_ncols);
        exit(-1);
        //return MORSE_ERR_ILLEGAL_VALUE; 
    }
    __svd(
            _M,
            maxrk,
            _CU,  ld_CU,
            _CV,  ld_CV, CU_ncols,
            newU,   ld_newU,
            newV,   ld_newV,  &new_UVrk ,
            rank, acc,
            svd_rA, svd_rB, svd_T, svd_sigma, svd_superb
         );
    if(use_scratch == 0) {
        free(svd_rA);
        free(svd_superb);
        free(svd_sigma);
    }
    if(use_scratch == 0 && use_trmm == 0){
        free(svd_rB);
    //    free(svd_T);
    }

    /** parameters that will be passed to HiCMA_HCORE_dgemm_ormqr */
    *_p_work_new     =  work         ;
    *_p__CU          =  _CU          ;
    *_p__CV          =  _CV          ;
    *_p_CU_ncols     =  CU_ncols     ;  
    *_p_new_UVrk     =  new_UVrk     ;
    *_p_newU         =  newU         ;
    *_p_ld_newU      =  ld_newU      ;
    *_p_qrtauA       =  qrtauA       ;
    *_p_CV_ncols     =  CV_ncols     ;
    *_p_newV         =  newV         ;
    *_p_ld_newV      =  ld_newV      ;
    *_p_qrtauB       =  qrtauB       ;
    *_p_use_CUV_clone=  use_CUV_clone;
    *_p_CUclone      =  CUclone      ;
    *_p_ld_CUclone   =  ld_CUclone   ;
    *_p__CU_save     =  _CU_save     ;
    *_p_CVclone      =  CVclone      ;
    *_p_ld_CVclone   =  ld_CVclone   ;
    *_p__CV_save     =  _CV_save     ;
}
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
)
{
    if(gemm_print_index){
		printf("%d:%s work:%p ", __LINE__, __func__, work);
		printf("M:%d N:%d  LDA:%d LDB:%d LDC:%d rk:%d maxrk:%d acc:%e a:%e b:%e\n",
				M, N, LDA, LDB, LDC, rk, maxrk, acc, alpha, beta);
	}

    int new_Crk = 0;
    /*
     * NOTES:
     * assumptions on matrix dimensions are marked as //ASSUMPTION
     */
    /*printf("%d %d|%g->%d. %g %g, %g %g\n",  */
            /*__LINE__, __COUNTER__,*/
            /*Crk[0], new_Crk, CU[0], CU[1], CV[0], CV[1]);*/
    //hcore_dgemm(Aij, ik, Ajk, -1, rank, acc);
    int _Ark = (int)(Ark[0]);
    int _Brk = (int)(Brk[0]);
    int _Crk = (int)(Crk[0]);
    /*if(_Ark == 0 || _Brk == 0 || _Crk == 0){
        fprintf(stderr, "%s %d: _Ark=%d _Brk=%d _Crk=%d. These rank values should not be zero.\n", __FILE__, __LINE__, _Ark, _Brk, _Crk);
        exit(-1);
        //return MORSE_ERR_ILLEGAL_VALUE; 
    }*/

    int _M = M; int _N = N;
    double* _CU = CU; int ld_CU = LDC;
    double* _CV = CV; int ld_CV = LDC; int* pnew_Crk = &new_Crk;
    double* _AU = AU; int ld_AU = LDA;
    double* _AV = AV; int ld_AV = LDA;
    double* _BU = BU; int ld_BU = LDB;
    double* _BV = BV; int ld_BV = LDB;
    int rank = rk;


    //printf("%s %s %d maxrank=%d\n", __FILE__, __func__, __LINE__, maxrk);
    char chall = 'A';
    double* CUclone = NULL;
    size_t CUclone_nelm =  _M * 2 * maxrk;
    int ld_CUclone = _M;
    int use_CUV_clone = 1;
    if(use_CUV_clone == 1) {
        if(use_scratch){
            CUclone = work;
            work += CUclone_nelm;
        } else {
            CUclone = malloc(CUclone_nelm * sizeof(double));
        }
        LAPACK_dlacpy(&chall,
                &_M, &_Crk,
                _CU, &ld_CU,
                CUclone, &ld_CUclone);
    }

    double* CVclone = NULL;
    size_t CVclone_nelm = _M * 2 * maxrk;
    int ld_CVclone = _M;
    if(use_CUV_clone == 1) {
        if(use_scratch){
            CVclone = work;
            work += CVclone_nelm;
        } else {
            CVclone = malloc(CVclone_nelm * sizeof(double));
        }
        LAPACK_dlacpy(&chall,
                &_M, &_Crk,
                _CV, &ld_CV,
                CVclone, &ld_CVclone);
    }
    double* _CU_save = _CU;
    double* _CV_save = _CV;

    if(use_CUV_clone == 1) {
        _CU = CUclone;
        _CV = CVclone;
    }
    int nb    = _M;//ASSUMPTION
    double* qrtauA = NULL;
    size_t qrtauA_nelm = nb;
    if(use_scratch){
        qrtauA = work;
    } else {
        qrtauA = malloc(qrtauA_nelm * sizeof(double));
    }
    assert(qrtauA != NULL);
    double* qrtauB = NULL;
    size_t qrtauB_nelm = nb;
    if(use_scratch){
        qrtauB = work + qrtauA_nelm;
    } else {
        qrtauB = malloc(qrtauB_nelm * sizeof(double));
    }
    assert(qrtauB != NULL);



    int CU_ncols = 0;
    //Warning:In this function, there are assumptions
    //on leading dimensions and number of rows/cols
    __qra(_M, maxrk, _CU, ld_CU, _Crk, &CU_ncols, _AU, ld_AU, _Ark, alpha, beta, qrtauA);
    //output: _CU contains [_CU _AU]. use _Crk+_Ark as number of cols
    assert(CU_ncols == (_Crk + _Ark));

    int CV_ncols = 0;
    //Warning:In this function, there are assumptions
    //on leading dimensions and number of rows/cols
    double* qrb_aubut = NULL;
    size_t qrb_aubut_nelm = maxrk * maxrk;
    if(use_scratch){
       qrb_aubut = work + qrtauA_nelm + qrtauB_nelm;
    } else {
       qrb_aubut = malloc(qrb_aubut_nelm * sizeof(double));
    }
    __qrb(_M, maxrk, _CV, ld_CV, _Crk, &CV_ncols,
            _AV, ld_AV, _Ark, _BU, ld_BU, _BV, ld_BV, _Brk, qrtauB, qrb_aubut);
    if(CU_ncols == 0 || CV_ncols == 0){
        fprintf(stderr, "%s %d: CU_ncols=%d CV_ncols=%d. These values should not be zero.\n", __FILE__, __LINE__, CU_ncols, CV_ncols);
        exit(-1);
        //return MORSE_ERR_ILLEGAL_VALUE; 
    }
    if(use_scratch == 0) {
        free(qrb_aubut);
    }
    assert(CU_ncols == CV_ncols);


    int ld_newU = nb; //ASSUMPTION
    int ld_newV = 2*maxrk; //ASSUMPTION
    int new_UVrk;
    double* newU = NULL;
    size_t newU_nelm = ld_newU * 2*maxrk;
    if(use_scratch){
        newU = work + qrtauA_nelm + qrtauB_nelm + qrb_aubut_nelm;
    } else {
        newU = malloc(newU_nelm * sizeof(double));
    }
    double* newV = NULL;
    size_t newV_nelm = nb * ld_newV;

    if(use_scratch){
        newV = work + qrtauA_nelm + qrtauB_nelm + qrb_aubut_nelm + newU_nelm;
    } else {
        newV = malloc(newV_nelm *  sizeof(double));
    }
    assert(newU != NULL);
    assert(newV != NULL);
    double *svd_rA = NULL;
    int svd_rA_nrows  = _M < CU_ncols ? _M : CU_ncols;
    int svd_rA_ncols  = CU_ncols;    //ASSUMPTION
    size_t svd_rA_nelm;
    if(use_trmm == 1){ // allocate more because rA will be output of trmm(rA*rB^T)
        svd_rA_nelm = svd_rA_nrows * svd_rA_nrows;
    } else {
        svd_rA_nelm = svd_rA_nrows * svd_rA_ncols;
    }
    if(use_scratch == 1) {
        svd_rA = work + qrtauA_nelm + qrtauB_nelm + qrb_aubut_nelm + newU_nelm + newV_nelm;
    } else {
        svd_rA = malloc(svd_rA_nelm * sizeof(double));
    }
    double *svd_rB = NULL;
    int svd_rB_nrows = _M < CV_ncols ? _M : CV_ncols;
    int svd_rB_ncols = CV_ncols;
    size_t svd_rB_nelm;
    if(use_trmm == 0) {
        svd_rB_nelm = svd_rB_nrows * svd_rB_ncols;
        if(use_scratch == 1) {
            svd_rB = work + qrtauA_nelm + qrtauB_nelm + qrb_aubut_nelm + newU_nelm + newV_nelm + svd_rA_nelm;
        } else {
            svd_rB = malloc(svd_rB_nelm * sizeof(double));
        }
    } else {
        svd_rB_nelm = 0; //This variable will be used for calculating offset for work array
    }
    double *svd_T = NULL;
    int svd_T_nrows = svd_rA_ncols;
    int svd_T_ncols = svd_rA_ncols;
    size_t svd_T_nelm;
    if(use_trmm == 0) {
        svd_T_nelm = svd_T_nrows * svd_T_ncols;
        if(use_scratch == 1) {
            svd_T = work + qrtauA_nelm + qrtauB_nelm + qrb_aubut_nelm + newU_nelm + newV_nelm + svd_rA_nelm + svd_rB_nelm;
        } else {
            svd_T = malloc(svd_T_nelm * sizeof(double));
        }
    } else {
        svd_T_nelm = 0; //This variable will be used for calculating offset for work array
    }

    double *svd_sigma  = NULL;
    double *svd_superb = NULL;
    size_t svd_sigma_nelm  = svd_T_nrows;
    size_t svd_superb_nelm = svd_T_nrows;
    if(use_scratch == 1){
        svd_sigma  = work + qrtauA_nelm + qrtauB_nelm + qrb_aubut_nelm + newU_nelm + newV_nelm + svd_rA_nelm + svd_rB_nelm + svd_T_nelm;
        svd_superb = work + qrtauA_nelm + qrtauB_nelm + qrb_aubut_nelm + newU_nelm + newV_nelm + svd_rA_nelm + svd_rB_nelm + svd_T_nelm + svd_sigma_nelm;
    } else {
        svd_sigma  = malloc(svd_sigma_nelm  * sizeof(double));
        svd_superb = malloc(svd_superb_nelm * sizeof(double));
    }
    if(ld_newV < CU_ncols){ /* ld_newV=2*maxrank CU_ncols=Crk+Ark */
        fprintf(stderr, "%s %d: Increase computation maxrank. %d is not enough ld_newV:%d CU_ncols:%d\n", __FILE__, __LINE__, maxrk, ld_newV, CU_ncols);
        exit(-1);
        //return MORSE_ERR_ILLEGAL_VALUE; 
    }
    __svd(
            _M,
            maxrk,
            _CU,  ld_CU,
            _CV,  ld_CV, CU_ncols,
            newU,   ld_newU,
            newV,   ld_newV,  &new_UVrk ,
            rank, acc,
            svd_rA, svd_rB, svd_T, svd_sigma, svd_superb
         );
    if(use_scratch == 0) {
        free(svd_rA);
        free(svd_superb);
        free(svd_sigma);
    }
    if(use_scratch == 0 && use_trmm == 0){
        free(svd_rB);
    //    free(svd_T);
    }

    /** parameters that will be passed to HiCMA_HCORE_dgemm_ormqr */
    *_p_work_new     =  work         ;
    *_p__CU          =  _CU          ;
    *_p__CV          =  _CV          ;
    *_p_CU_ncols     =  CU_ncols     ;  
    *_p_new_UVrk     =  new_UVrk     ;
    *_p_newU         =  newU         ;
    *_p_ld_newU      =  ld_newU      ;
    *_p_qrtauA       =  qrtauA       ;
    *_p_CV_ncols     =  CV_ncols     ;
    *_p_newV         =  newV         ;
    *_p_ld_newV      =  ld_newV      ;
    *_p_qrtauB       =  qrtauB       ;
    *_p_use_CUV_clone=  use_CUV_clone;
    *_p_CUclone      =  CUclone      ;
    *_p_ld_CUclone   =  ld_CUclone   ;
    *_p__CU_save     =  _CU_save     ;
    *_p_CVclone      =  CVclone      ;
    *_p_ld_CVclone   =  ld_CVclone   ;
    *_p__CV_save     =  _CV_save     ;
}

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
)
{
    char chall = 'A';
    int new_Crk = new_UVrk;
    int _Crk = (int)(Crk[0]);
    int _M = M; int _N = N;
    int ld_CU = LDC;
    int ld_CV = LDC;
    int ncols_qA = CU_ncols;
    if(new_UVrk > storage_maxrank) {
        fprintf(stderr, "%s %d: Increase storage maxrank %d. Observed rank for U and V is %d\n", __FILE__, __LINE__, storage_maxrank, new_UVrk);
        exit(-1);
    }
    __newu(_M,
            ncols_qA,
            _CU, ld_CU, _Crk,
            newU, ld_newU, new_UVrk,
            qrtauA
            );
    int ncols_qB = CV_ncols;
    __newv(_M,
            ncols_qB,
            _CV, ld_CV, _Crk,
            newV, ld_newV, new_UVrk,
            qrtauB
            );
    if(use_CUV_clone == 1) {
        LAPACK_dlacpy(&chall,
                &_M, &new_UVrk,
                CUclone, &ld_CUclone,
                _CU_save, &ld_CU
                );

        LAPACK_dlacpy(&chall,
                &_M, &new_UVrk,
                CVclone, &ld_CVclone,
                _CV_save, &ld_CV
                );
        if(use_scratch == 0) {
            free(CUclone);
            free(CVclone);
        }
    }
    if(use_scratch == 0) {
        free(qrtauA);
        free(qrtauB);
        free(newU);
        free(newV);
    }
    int old_Crk = Crk[0];
    if(gemm_print_index){
        printf("Crk[0]:%d %g RANK CHANGE: %d->%d\n",
                Crk[0], Crk[0], old_Crk, new_Crk);
    }
    Crk[0] = new_Crk;
    if(gemm_print_index){
        int casted_Crk = (int)(Crk[0]);
        printf("casted_Crk:%d  Crk[0]:%d %g RANK CHANGE: %d->%d\n",
                casted_Crk, Crk[0], Crk[0], old_Crk, new_Crk);
    }
}

double
HiCMA_check_norm(int _n_, int _ld_, double* _A_, double* _B_, double* p_elem_work){
    /* Copy _B_ into p_work */
    cblas_dcopy(_n_*_n_, _B_, 1, p_elem_work, 1);

    /*((double*)p_elem_work)[0] = 4.0;*/
    /*printf("%g\n", ((double*)p_elem_work)[0]);*/

    /* Check _A_ and _B_ are equal */
    /* ~/parsec/build/dplasma/cores/core_dgeadd.c */
    /* void cblas_daxpy (const MKL_IN_A_ n, const double a, const double *x, const MKL_IN_A_ incx, double *y, const MKL_IN_A_ incy); */
    /* y := a*x + y */
    /* _B_ = -1 * _A_ + _B_ */
    cblas_daxpy(_n_*_n_, (double)-1., _A_, 1, p_elem_work, 1);
    double norm = 0.0;
    /*CORE_dlansy(PlasmaInfNorm, uplo, _n_, p_elem_work, _ld_, p_elem_norm, &norm); */
    //norm = LAPACKE_dlansy_work(LAPACK_COL_MAJOR, 'I', 'L', _n_, p_elem_work, _ld_, p_elem_norm);
    norm = LAPACKE_dlansy(LAPACK_COL_MAJOR, 'I', 'L', _n_, p_elem_work, _ld_); // FIXME CORE_d... function must be used 
    if(norm > HiCMA_HCORE_ERROR_THRESHOLD_DENSE) {
        printf("norm %g is larger than HiCMA_HCORE_ERROR_THRESHOLD_DENSE %g. Exiting...\n", norm, HiCMA_HCORE_ERROR_THRESHOLD_DENSE);
    //    exit(-1);
    }
    return norm;
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

/* Used for MPI_Ggatherv */
void HICMA_get_stat2(int *G, size_t m, int maxrank, HICMA_stat_t *stat)
{
    int min = maxrank, max = 0, sum = 0;
    long int num = m;

    for(int i = 0; i < m; i++){
        if( G[i] <= maxrank && G[i] >= 0 ){
            sum += G[i];
            if( G[i] < min )
                min = G[i];
            if( G[i] > max )
                max = G[i];
        } else {
            num -= 1;
        } 
    }
    printf("number of tiles:%lu\n", num);
    fflush(stdout);

    stat->min = min;
    stat->max = max;
    stat->avg = ((double)sum)/num;
}

void HICMA_print_stat(HICMA_stat_t stat)
{
    printf("avg:%g min:%d max:%d\n", stat.avg, stat.min, stat.max);
}
