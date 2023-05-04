/**
 * @copyright (c) 2021     King Abdullah University of Science and Technology (KAUST).
 * @copyright (c) 2021     The Universiy of Tennessee and The Universiy of Tennessee Research Foundation.
 *                         All rights reserved.
 * @version 0.1.0
 * @date 2021-01-24
 **/

#include "hicma_parsec.h"

/* Problem types */
char* str_problem[8]={"randtlr", "ed-2d-sin", "st-2d-sqexp", "st-3d-sqexp", "st-3d-exp", "ed-3d-sin",  "md-3d-virus", "md-3d-cube"};

/* Used for gather operation count */
unsigned long *op_band, *op_offband, *op_path, *op_offpath;

/* Kernel of generating matrix */
static int starsh_generate_map_operator(parsec_execution_stream_t *es,
        const parsec_tiled_matrix_dc_t *descA, void *_A, int uplo,
        int m, int n, void *op_data)
{
    int tempmm, tempnn, ldam;
    double *A = _A;
    starsh_params_t *params = op_data;
    (void)es;
    (void)uplo;

    tempmm = ((m)==((descA->mt)-1)) ? ((descA->m)-(m*(descA->mb))) : (descA->mb);
    tempnn = ((n)==((descA->nt)-1)) ? ((descA->n)-(n*(descA->nb))) : (descA->nb);
    ldam   = BLKLDD( descA, m );

    params->kernel(tempmm, tempnn, params->index + m*descA->mb,
            params->index + n*descA->nb, params->data, params->data, A, ldam);

    return 0;
}

/* Generate metrix */
static int starsh_generate_map(parsec_context_t *parsec, int uplo,
        parsec_tiled_matrix_dc_t *A, starsh_params_t *params)
{
    parsec_taskpool_t *pool = NULL;

    /* Check input arguments */
    if ((uplo != dplasmaLower) &&
        (uplo != dplasmaUpper) &&
        (uplo != dplasmaUpperLower))
    {
        dplasma_error("starsh_generate_map", "illegal value of type");
        return -3;
    }

    starsh_params_t *copy_params = malloc(sizeof(starsh_params_t));
    *copy_params = *params;
    pool = parsec_apply_New(uplo, A, starsh_generate_map_operator, copy_params);

    if ( pool != NULL ) {
        parsec_context_add_taskpool(parsec, pool);
        parsec_context_start(parsec);
        parsec_context_wait(parsec);
        parsec_apply_Destruct(pool);
    }
    return 0;
}

/* Rank statistics */
void tlr_rank_stat(char* strid, int *rank_array, int band_size,  parsec_context_t* parsec, two_dim_block_cyclic_t* pdcAr, int* pminrk, int* pmaxrk, double* pavgrk, int maxrank){
    two_dim_block_cyclic_t dcAr = *pdcAr;
    int NT = dcAr.super.lm;
    int rank = dcAr.super.super.myrank;
    int root = dcAr.super.super.rank_of(&dcAr.super.super, 0, 0);

    /* Timer start */
    SYNC_TIME_START();

    /* Gather from dcAr to G (global) */
    parsec_rank_gather(parsec, (parsec_tiled_matrix_dc_t*)&dcAr, rank_array, band_size);

    /* Calculate min/avg/max */
    HICMA_stat_t rankstat;
    HICMA_get_stat2(rank_array, dcAr.super.lm, band_size, &rankstat);

    *pminrk = rankstat.min;
    *pmaxrk = rankstat.max;
    *pavgrk = rankstat.avg;

    /* Timer end */
    SYNC_TIME_PRINT(rank, ("Gather rank" "\tband_size= %d  NT= %4d avg= %lf min= %d max= %d\n",
			    band_size, NT, *pavgrk, *pminrk, *pmaxrk));

#if PRINT_RANK 
    SYNC_TIME_START();
    if(rank == 0){ // Takes too much time for big matrices
	    printf("%s %d %d\n", strid, NT, NT);
	    int i, j;
	    for(i = 0; i < dcAr.super.lm; i++){
		    for(j = 0; j < dcAr.super.ln; j++){//ASSUMPTION MT==NT
			    if( i - j >= band_size ) 
				    printf("%3d ", rank_array[j*dcAr.super.lm+i]);
			    else
				    printf("%3d ", -1); //%-3d
		    }
		    printf("\n");
	    }
    }
    SYNC_TIME_PRINT(rank, ("Print ranks of tiles\n"));
#endif
}

/* Check difference of ||L*L'-A|| */
int check_dpotrf2( parsec_context_t *parsec, int verbose,
		int uplo,
		parsec_tiled_matrix_dc_t *A,
		parsec_tiled_matrix_dc_t *A0, double threshold )
{
    two_dim_block_cyclic_t *twodA = (two_dim_block_cyclic_t *)A0;
    two_dim_block_cyclic_t LLt;
    int info_factorization;
    double Rnorm = 0.0;
    double Anorm = 0.0;
    double result = 0.0;
    int N = A->n;
    double eps = LAPACKE_dlamch_work('e');
    int side;

    two_dim_block_cyclic_init(&LLt, matrix_RealDouble, matrix_Tile,
                              A->super.nodes, twodA->grid.rank,
                              A->mb, A->nb, N, N, 0, 0, N, N,
                              twodA->grid.krows, twodA->grid.kcols,
                              twodA->grid.rows);

    LLt.mat = parsec_data_allocate((size_t)LLt.super.nb_local_tiles *
                                  (size_t)LLt.super.bsiz *
                                  (size_t)parsec_datadist_getsizeoftype(LLt.super.mtype));

    dplasma_dlaset( parsec, dplasmaUpperLower, 0., 0.,(parsec_tiled_matrix_dc_t *)&LLt );
    dplasma_dlacpy( parsec, uplo, A, (parsec_tiled_matrix_dc_t *)&LLt );

    /* Compute LL' or U'U  */
    side = (uplo == dplasmaUpper ) ? dplasmaLeft : dplasmaRight;
    dplasma_dtrmm( parsec, side, uplo, dplasmaTrans, dplasmaNonUnit, 1.0,
                   A, (parsec_tiled_matrix_dc_t*)&LLt);

    /* compute LL' - A or U'U - A */
    dplasma_dtradd( parsec, uplo, dplasmaNoTrans,
                    -1.0, A0, 1., (parsec_tiled_matrix_dc_t*)&LLt);

    Anorm = dplasma_dlansy(parsec, dplasmaFrobeniusNorm, uplo, A0);
    Rnorm = dplasma_dlansy(parsec, dplasmaFrobeniusNorm, uplo,
                           (parsec_tiled_matrix_dc_t*)&LLt);

    //result = Rnorm / ( Anorm * N * eps ) ;
    result = Rnorm / ( Anorm ) ;

    if ( verbose ) {
        printf("============\n");
        printf("Checking the Cholesky factorization. Threshold is " BLU "%.2e" RESET " \n", threshold);
        printf( "-- ||A||_oo = %e, ||L'L-A||_oo = %e\n", Anorm, Rnorm );
        printf("-- ||L'L-A||_oo/(||A||_oo) = " RED "%e" RESET "\n", result);
    }

    if ( isnan(Rnorm)  || isinf(Rnorm)  ||
         isnan(result) || isinf(result) ||
         //(result > 60.0) )
         (result > threshold) )
    {
        if( verbose ) printf(RED "-- Factorization is suspicious ! \n" RESET);
        info_factorization = 1;
    }
    else
    {
        if( verbose ) printf(GRN "-- Factorization is CORRECT ! \n" RESET);
        info_factorization = 0;
    }

    parsec_data_free(LLt.mat); LLt.mat = NULL;
    parsec_tiled_matrix_dc_destroy( (parsec_tiled_matrix_dc_t*)&LLt);

    return info_factorization;
}

/* Check solution L difference from the dense counterpart */ 
int check_diff( parsec_context_t *parsec, int verbose,
                  int uplo,
                  parsec_tiled_matrix_dc_t *A,
                  parsec_tiled_matrix_dc_t *A0, double threshold )
{
    two_dim_block_cyclic_t *twodA = (two_dim_block_cyclic_t *)A0;
    two_dim_block_cyclic_t LLt;
    int info_factorization;
    double Rnorm = 0.0;
    double A0norm = 0.0;
    double Anorm = 0.0;
    double result = 0.0;
    int N = A->n;
    double eps = LAPACKE_dlamch_work('e');
    int side;

    Anorm = dplasma_dlansy(parsec, dplasmaFrobeniusNorm, uplo, A);
    A0norm = dplasma_dlansy(parsec, dplasmaFrobeniusNorm, uplo, A0);

    /* compute A = A - A0 */
    dplasma_dtradd( parsec, uplo, dplasmaNoTrans,
                    -1.0, A0, 1., A);

    Rnorm = dplasma_dlansy(parsec, dplasmaFrobeniusNorm, uplo, A);

    if ( verbose ) {
        printf("============\n");
        printf("Checking the equality of lower/upper triangular parts of two matrices. Threshold is %.2e\n", threshold);
        printf( "-- dpotrf: ||L0||_oo = %e, HiCMA: ||L||_oo = %e, ||L-L0||_oo = %e, ||L-L0||_oo/||L||_oo = %e \n", A0norm, Anorm, Rnorm, Rnorm/Anorm );
    }

    if ( isnan(Rnorm)  || isinf(Rnorm)  ||
         isnan(result) || isinf(result) ||
         //(result > 60.0)
         (Rnorm/Anorm > threshold)
         )
    {
        if( verbose ) printf("-- DIFFERENT matrices ! \n");
        info_factorization = 1;
    }
    else
    {
        if( verbose ) printf("-- SAME matrices ! \n");
        info_factorization = 0;
    }

    return info_factorization;
}

#if defined(PARSEC_HAVE_LCI)
static void lci_sum_op_ul(void *dst, const void *src, size_t count)
{
    unsigned long *d = dst;
    const unsigned long *s = src;
    size_t c = count / sizeof(unsigned long);
    for (size_t i = 0; i < c; i++)
        d[i] += s[i];
}

static void lci_max_op_ul(void *dst, const void *src, size_t count)
{
    unsigned long *d = dst;
    const unsigned long *s = src;
    size_t c = count / sizeof(unsigned long);
    for (size_t i = 0; i < c; i++)
        if (s[i] > d[i])
            d[i] = s[i];
}

static void lci_sum_op_d(void *dst, const void *src, size_t count)
{
    double *d = dst;
    const double *s = src;
    size_t c = count / sizeof(double);
    for (size_t i = 0; i < c; i++)
        d[i] += s[i];
}
#endif

int main(int argc, char ** argv)
{
    parsec_context_t* parsec;
    int iparam[IPARAM_SIZEOF] = {0};
    double dparam[DPARAM_SIZEOF];
    int uplo = dplasmaLower;
    int info = 0;
    int ret = 0;
    double time_opt_band = 0.0, time_regenerate = 0.0;

    /* Initialize PaRSEC */
    parsec = setup_parsec(argc, argv, iparam, dparam);
    PASTE_CODE_IPARAM_LOCALS(iparam);
    PASTE_CODE_DPARAM_LOCALS(dparam);
    PASTE_CODE_FLOPS(FLOPS_DPOTRF, ((double)N));

    if( 0 == rank && verbose ) { printf("%s is starting\n", argv[0]); fflush(stdout); }
    if( 0 == rank && verbose ) {
        printf("M:%d N:%d NB:%d NB:%d HNB:%d HNB:%d\n", N, N, NB, NB, HNB, HNB);
        printf("nodes:%d P:%d Q:%d cores:%d\n", nodes, P, Q, cores);
        printf("kind_of_problem:%d %s\n", kind_of_problem, str_problem[kind_of_problem]);
        printf("tol:%.1e add_diag:%g fixed_rk:%d wave_k:%g\n", tol, add_diag, fixedrk, wave_k);
        printf("send_full_tile:%d lookahead: %d band_size:%d max_rank:%d gen:%d comp:%d\n", send_full_tile, lookahead, band_size, maxrank, genmaxrank, compmaxrank);
    }

    /* Band distribution */
    int P_BAND = 1;
    int NB_BAND = NB;    // the rank on band
    int N_BAND = NB_BAND * NT;

    /* Off band */
    int NB_UV = NB;      // NB the same as band, because subtile_desc_create will pass desc->nb
    int N_UV = NB_UV * NT;

    /* Indicator for free memory */
    int indicator_band = 0;
    int indicator_offband = 1;

    /* init STARS-H data to fill matrix */
    void *data;
    STARSH_kernel *kernel;
    
    /* Default parameters for statistics */
    double beta = 0.1;
    double nu = 0.5;     //in matern, nu=0.5 exp (half smooth), nu=inf sqexp (inifinetly smooth)
    double noise = 1.e-1;
    double sigma = 1.0;
    /* Placement template for particles */
    enum STARSH_PARTICLES_PLACEMENT place = STARSH_PARTICLES_UNIFORM;
        
    srand(0);

    /* Synthetic matrix with equal ranks */
    if(kind_of_problem == 0)
    {
        kernel = starsh_randtlr_block_kernel;
        double decay = 0.5;
        info = starsh_randtlr_generate((STARSH_randtlr **)&data, N, NB, decay,
                add_diag);
    }
    /* Semi-synthetic matrix with non-equal ranks
     * electrodynamics application from STARS-H
     */
    else if(kind_of_problem == 1)
    {
        /* Dimensionality of the problem (do not change it) */
        int ndim = 2; // 2 is used in europar   3;
        /* Kernel must be in accordance with problem dimensionality */
        kernel = starsh_eddata_block_sin_kernel_2d;
        /* Generating data */
        info = starsh_eddata_generate((STARSH_eddata **)&data, N, ndim, wave_k,
                add_diag, place);
    } 
    /* statistics-2d-sqexp */
    else if(kind_of_problem == 2) {
        int ndim = 2;
        kernel = starsh_ssdata_block_sqrexp_kernel_2d; 
        info = starsh_ssdata_generate((STARSH_ssdata **)&data, N, ndim,
                beta, nu, noise,
                place, sigma);
    }
    /* statistics-3d-sqexp */
    else if(kind_of_problem == 3) {
        int ndim = 3;
        kernel = starsh_ssdata_block_sqrexp_kernel_3d; 
        info = starsh_ssdata_generate((STARSH_ssdata **)&data, N, ndim,
                beta, nu, noise,
                place, sigma);
    }
    /* statistics-3d-exp */
    else if(kind_of_problem == 4) {
        int ndim = 3;
        kernel = starsh_ssdata_block_exp_kernel_3d; 
        info = starsh_ssdata_generate((STARSH_ssdata **)&data, N, ndim,
                beta, nu, noise,
                place, sigma);
    }
    /* electrodynamics-3d-sin */
    else if(kind_of_problem == 5) {
        int ndim = 3;
        kernel = starsh_eddata_block_sin_kernel_3d;
        info = starsh_eddata_generate((STARSH_eddata **)&data, N, ndim, wave_k,
                add_diag, place);
    }
    else if(kind_of_problem == 6){
        int ndim = 3;
        printf("\nRBF, %s, %d, %d, %d, %d, %f, %f\n", mesh_file, N, rbf_kernel, numobj, order, radius, density);
        kernel = starsh_generate_3d_virus;
	starsh_generate_3d_rbf_mesh_coordinates_virus((STARSH_mddata **)&data, 
                mesh_file, N, ndim, rbf_kernel, numobj, 1, add_diag, 
                radius, density, order);

    }
    else if(kind_of_problem == 7){
        int ndim = 3;
        kernel = starsh_generate_3d_cube;
        starsh_generate_3d_rbf_mesh_coordinates_cube((STARSH_mddata **)&data, 
                N, ndim, rbf_kernel, 1, add_diag, radius, order);

    }
    else
    {
        printf("Wrong value of \"kind_of_problem\" parameter\n");
        return -1;
    }

    if(info != 0)
    {
        printf("Problem was NOT generated (wrong parameters)\n");
        return info;
    }

    STARSH_int *index = malloc(N*sizeof(STARSH_int));
    STARSH_int i; 
    for(i = 0; i < N; ++i)
        index[i] = i;
    starsh_params_t params = {data, kernel, index};

    /* Make sure lookahead > 0 */
    assert(lookahead >= -1);

    /* Make sure band >= 0 */
    assert(band_size >= 1);

    /* Make sure NB_BAND >= NB_UV, meaning NB >= 2 * maxrank */
    assert(NB_BAND >= NB_UV);

    /* If auto select band_size, band_size set to 1 at the beginning */
    if( auto_band ) {
        if( rank == 0 && verbose) { 
            printf(YEL "\n%d: Auto-tuning band_size, set band_size = 1 at the beginning\n" RESET, __LINE__);
            fflush(stdout);
        }
        band_size = 1;
    }

    /* dcA data descriptor */
    sym_two_dim_block_cyclic_band_t dcA;
    sym_two_dim_block_cyclic_init(&dcA.off_band, matrix_RealDouble,
                                  nodes, rank, NB, NB_UV, N, N_UV, 0, 0,
                                  N, N_UV, P, uplo);

    /* Init band */
    two_dim_block_cyclic_init(&dcA.band, matrix_RealDouble, matrix_Tile,
                              nodes, rank, NB, NB_BAND, NB*band_size, N_BAND,
                              0, 0, NB*band_size, N_BAND, 1, 1, P_BAND);
#if BAND_MEMORY_CONTIGUOUS
    /* Allocate memory on band */
    dcA.band.mat = parsec_data_allocate((size_t)dcA.band.super.nb_local_tiles *
                                   (size_t)dcA.band.super.bsiz *
                                   (size_t)parsec_datadist_getsizeoftype(dcA.band.super.mtype));
#endif

    /* Init two_dim_block_cyclic_band_t structure */
    hicma_parsec_sym_two_dim_block_cyclic_band_init( &dcA, nodes, rank, band_size );
    parsec_data_collection_set_key((parsec_data_collection_t*)&dcA, "dcA_off_band");
    parsec_data_collection_set_key(&dcA.band.super.super, "dcA_band");


    /* dcAr contains rank of each tile. It is NT by NT matrix for now. */
    two_dim_block_cyclic_t dcAr;
    two_dim_block_cyclic_init(&dcAr, matrix_Integer, matrix_Tile,
                              nodes, rank, 1, 1, NT, NT, 0, 0,
                              NT, NT, 1, 1, P);
    dcAr.mat = parsec_data_allocate((size_t)dcAr.super.nb_local_tiles *
                                   (size_t)dcAr.super.bsiz *
                                   (size_t)parsec_datadist_getsizeoftype(dcAr.super.mtype));
    parsec_data_collection_set_key((parsec_data_collection_t*)&dcAr, "dcAr");

    /* Initialize dcAr to negtive, used to gather rank using MPI_Gather / MPI_Gatherv */
    for(int i = 0; i < dcAr.super.nb_local_tiles; i++) 
        ((int *)dcAr.mat)[i] = -1; 

    if( 0 == rank && verbose ) { printf("%d: Dense diagonal, U, V, rank Matrices are allocated\n", __LINE__);fflush(stdout);}
    if( 0 == rank && verbose ) { printf("%d: STARSH will generate problem\n", __LINE__);fflush(stdout);}

    /* Timer start */
    SYNC_TIME_START();

    /* matrix generation */
    STARSH_gen(parsec, uplo, (parsec_tiled_matrix_dc_t*)&dcA,
                       (parsec_tiled_matrix_dc_t*)&dcAr,
                       data, kernel, index, tol, genmaxrank,
		       band_size, send_full_tile, &info);

    /* Timer end */
    SYNC_TIME_PRINT(rank, ("STARSH_gen" "\tPxQ= %3d %-3d NB= %4d N= %7d\n",
                           P, Q, NB, N)); 

    double time_starsh = sync_time_elapsed;

    if(info != 0) {
        printf("Error in low rank matrix generation, info:%d\n", info);
        fflush(stdout);
        cleanup_parsec(parsec, iparam);
        return -1;
    }

#if DEBUG_INFO
    /* Timer start */
    SYNC_TIME_START();

    /* Check correctness of rank */
    parsec_rank_check(parsec, (parsec_tiled_matrix_dc_t*)&dcAr, band_size);

    /* Timer end */
    SYNC_TIME_PRINT(rank, ("Before HiCMA Check Rank" "\tPxQ= %3d %-3d NB= %4d N= %7d maxrank= %d "
                           "band_size= %d send_full= %d\n\n",
                           P, Q, NB, N, maxrank, band_size, send_full_tile));
#endif

    /* Gathering rank info */
    int imaxrk = -1, iminrk=-1;
    double iavgrk = -1.0;

    /* Used for gathering rank */
    int *rank_array = (int *)malloc(dcAr.super.lmt * dcAr.super.lnt * sizeof(int));

#if GATHER_RANK
    tlr_rank_stat("init_rank_tile", rank_array, band_size, parsec, &dcAr, &iminrk, &imaxrk, &iavgrk, maxrank);
#else
    /* Disable band_size auto-tuning, becasue the rank info is not gathered */
    auto_band = 0;
#endif

    int imaxrk_before = imaxrk;
    int iminrk_before = iminrk;
    double iavgrk_before = iavgrk;

    if( auto_band ) {
        /* Make sure band_size starts from 1 */
        assert( band_size == 1 );

        /* dcFake to control process rank working on*/ 
        two_dim_block_cyclic_t dcFake;
        two_dim_block_cyclic_init(&dcFake, matrix_Integer, matrix_Tile,
                                  nodes, rank, 1, 1, 1, nodes, 0, 0,
                                  1, nodes, 1, 1, 1);
        dcFake.mat = parsec_data_allocate((size_t)dcFake.super.nb_local_tiles *
                                          (size_t)dcFake.super.bsiz *
                                          (size_t)parsec_datadist_getsizeoftype(dcFake.super.mtype));
        parsec_data_collection_set_key((parsec_data_collection_t*)&dcFake, "dcFake");

        /* Timer start */
        SYNC_TIME_START();

        /* Find the best band_size */
        int band_size_opt = parsec_band_size_auto_tuning(parsec, (parsec_tiled_matrix_dc_t*)&dcAr,
                                                    (parsec_tiled_matrix_dc_t*)&dcFake,
                                                    rank_array, NB);

        /* Timer end */
        SYNC_TIME_PRINT(rank, ("band_size_auto_tuning" "\tPxQ= %3d %-3d NB= %4d N= %7d maxrank= %d "
                               "send_full= %d band_size_opt= %d\n\n",
                               P, Q, NB, N, maxrank, send_full_tile, band_size_opt));

	time_opt_band = sync_time_elapsed;

        /* Free memory */
        parsec_tiled_matrix_dc_destroy( (parsec_tiled_matrix_dc_t*)&dcFake );

        /* Timer start */
        SYNC_TIME_START();

        /* If band_size changed */
        if( band_size_opt != 1 ) {
            /* Free band memory */
            if( NULL == dcA.band.mat )
                parsec_band_free(parsec, (parsec_tiled_matrix_dc_t *)&dcA, band_size, indicator_band);
            else
                parsec_data_free(dcA.band.mat);
            parsec_tiled_matrix_dc_destroy( &dcA.band.super );

            /* Set band_size */
            band_size = band_size_opt; 

            /* Re-init band */
            two_dim_block_cyclic_init(&dcA.band, matrix_RealDouble, matrix_Tile,
                                      nodes, rank, NB, NB_BAND, NB*band_size, N_BAND,
                                      0, 0, NB*band_size, N_BAND, 1, 1, P_BAND);
#if BAND_MEMORY_CONTIGUOUS
            /* Re-allocate memory on band */
            dcA.band.mat = parsec_data_allocate((size_t)dcA.band.super.nb_local_tiles *
                                                (size_t)dcA.band.super.bsiz *
                                                (size_t)parsec_datadist_getsizeoftype(dcA.band.super.mtype));
#endif

	    /* Set band size */
            dcA.band_size = (unsigned int)band_size;

            /* band generation */
            parsec_band_regenerate(parsec, uplo, (parsec_tiled_matrix_dc_t *)&dcA,
			    data, kernel, index, band_size);
        }

        /* Timer end */
        SYNC_TIME_PRINT(rank, ("band_regeneration" "\tPxQ= %3d %-3d NB= %4d N= %7d maxrank= %d "
                               "send_full= %d band_size= %d\n\n",
                               P, Q, NB, N, maxrank, send_full_tile, band_size));

	time_regenerate = sync_time_elapsed;
    }


    /* Calculate memory needed for matrix allocation */
    /* Timer start */
    SYNC_TIME_START();

    double memory_per_node, memory_per_node_max;
    long long int *size_allocate = (long long int *)calloc(1, sizeof(long long int));
    long long int *size_allocate_max = (long long int *)calloc(1, sizeof(long long int));

    /* Calculate memory needed before factorization
     * memory_per_node: based on actual rank 
     * memory_per_node_max: based on maxrank
     */
    *size_allocate = (long long int)(NT - band_size) * (NT - band_size + 1) * NB * iavgrk + (long long int)band_size * NT * NB * NB;
    *size_allocate_max = (long long int)(NT - band_size) * (NT - band_size + 1) * NB * maxrank + (long long int)band_size * NT * NB * NB;

    memory_per_node = *size_allocate / (double)1024 / 1024 / 1024 * 8 / dcA.super.super.nodes;
    memory_per_node_max = *size_allocate_max / (double)1024 / 1024 / 1024 * 8 / dcA.super.super.nodes;

    /* Free memory */
    free(size_allocate);
    free(size_allocate_max);

    /* Timer end */
    SYNC_TIME_PRINT(rank, ("Memory_count" "\tPxQ= %3d %-3d NB= %4d N= %7d maxrank= %d "
                           "band_size= %d send_full= %d memory_per_node= %lf GB memory_per_node_max= %lf GB\n",
                           P, Q, NB, N, maxrank, band_size, send_full_tile, memory_per_node, memory_per_node_max));

    if( memory_per_node > THRESHOLD_MEMORY_PER_NODE ) {
            fprintf(stderr, "memory_for_matrix_allocation_per_node : %d < %lf Gbytes\n", THRESHOLD_MEMORY_PER_NODE, memory_per_node);
            return 0;
    }

    /* Gather rank info during Cholesky if PRINT_RANK*/
    /* dcRank data descriptor : init_rank, min_rank, max_rank, final_rank */
    sym_two_dim_block_cyclic_band_t dcRank;
    sym_two_dim_block_cyclic_init(&dcRank.off_band, matrix_Integer,
                                  nodes, rank, 1, RANK_MAP_BUFF,
                                  NT, RANK_MAP_BUFF*NT, 0, 0,
                                  NT, RANK_MAP_BUFF*NT, P, uplo);

    /* Init band */
    two_dim_block_cyclic_init(&dcRank.band, matrix_Integer, matrix_Tile,
                              nodes, rank, 1, RANK_MAP_BUFF,
                              band_size, RANK_MAP_BUFF*NT, 0, 0,
                              band_size, RANK_MAP_BUFF*NT, 1, 1, P_BAND);

    /* Init two_dim_block_cyclic_band_t structure */
    hicma_parsec_sym_two_dim_block_cyclic_band_init( &dcRank, nodes, rank, band_size );
    parsec_data_collection_set_key((parsec_data_collection_t*)&dcRank, "dcRank_super");
    parsec_data_collection_set_key(&dcRank.band.super.super, "dcRank_band");


    /* Used for checking results */
    sym_two_dim_block_cyclic_t dcAd;
    sym_two_dim_block_cyclic_t dcA0;

    if( !info && check)
    {
        /* dcAd */
        sym_two_dim_block_cyclic_init(&dcAd, matrix_RealDouble,
                                      nodes, rank, NB, NB, N, N, 0, 0,
                                      N, N, P, uplo);
        dcAd.mat = parsec_data_allocate((size_t)dcAd.super.nb_local_tiles *
                                        (size_t)dcAd.super.bsiz *
                                        (size_t)parsec_datadist_getsizeoftype(dcAd.super.mtype));
        parsec_data_collection_set_key((parsec_data_collection_t*)&dcAd, "dcAd");

        /* Generates original problem withOUT approximation in dcAd */
        starsh_generate_map(parsec, uplo, (parsec_tiled_matrix_dc_t *)&dcAd,
                &params);

        /* dcA0 */
        sym_two_dim_block_cyclic_init(&dcA0, matrix_RealDouble,
                                  nodes, rank, NB, NB, N, N, 0, 0,
                                  N, N, P, uplo);
        dcA0.mat = parsec_data_allocate((size_t)dcA0.super.nb_local_tiles *
                                        (size_t)dcA0.super.bsiz *
                                        (size_t)parsec_datadist_getsizeoftype(dcA0.super.mtype));
        parsec_data_collection_set_key((parsec_data_collection_t*)&dcA0, "dcA0");

        /* Timer start */
        SYNC_TIME_START();

        /* Uncompresses approximate matrix dcA into dcA0 */
        STARSH_check(parsec, uplo, (parsec_tiled_matrix_dc_t*)&dcA0,
                             (parsec_tiled_matrix_dc_t*)&dcA,
                             (parsec_tiled_matrix_dc_t*)&dcAr,
                             band_size, maxrank, &info);

        /* Timer end */
        SYNC_TIME_PRINT(rank, ("STARSH_check" "\tPxQ= %3d %-3d NB= %4d N= %7d\n",
                               P, Q, NB, N));

        if(rank == 0) {
            printf("Check matrix generation:\n");
        }
        double norm1 = dplasma_dlansy(parsec, dplasmaFrobeniusNorm, uplo,
                (parsec_tiled_matrix_dc_t *)&dcAd);
        double norm2 = dplasma_dlansy(parsec, dplasmaFrobeniusNorm, uplo,
                (parsec_tiled_matrix_dc_t *)&dcA0);
        if(rank == 0) {
            printf("Norms: map=%e lr_jdf=%e\n", norm1, norm2);
        }

        /* B = alpha * op(A) + beta * B*/
        /* dcA0 = -1 * dcAd + 1 * dcA0 */
        dplasma_dtradd(parsec, uplo, dplasmaNoTrans, -1.0,
                (parsec_tiled_matrix_dc_t *)&dcAd, 1.0,
                (parsec_tiled_matrix_dc_t *)&dcA0);

        double diff = dplasma_dlansy(parsec, dplasmaFrobeniusNorm, uplo,
                (parsec_tiled_matrix_dc_t *)&dcA0);

        if(rank == 0) {
            printf("Norms: map-lr_jdf=%e\n", diff);
        }

        if(rank == 0) {
            printf("Relative error of approximation: %e\n\n", diff/norm1);
        }

    }

    /* Used for operation count */
    unsigned long* tileopcounters = NULL; // count all operations performed on a tile
    int nelm_tileopcounters = NT*NT;
    tileopcounters = calloc(nelm_tileopcounters, sizeof(unsigned long)); //TODO free
    unsigned long*    opcounters  = NULL; // count all operations performed by a core
    int nelm_opcounters = nodes*cores;
    opcounters     = calloc(nelm_opcounters,     sizeof(unsigned long)); //TODO free

    /* band, offband, critical path, off critical path */
    op_band = (unsigned long *)calloc(cores, sizeof(unsigned long));
    op_offband = (unsigned long *)calloc(cores, sizeof(unsigned long));
    op_path = (unsigned long *)calloc(cores, sizeof(unsigned long));
    op_offpath = (unsigned long *)calloc(cores, sizeof(unsigned long));

#if DEBUG_INFO
    /* Timer start */
    SYNC_TIME_START();

    /* Check correctness of rank */
    parsec_rank_check(parsec, (parsec_tiled_matrix_dc_t*)&dcAr, band_size);

    /* Timer end */
    SYNC_TIME_PRINT(rank, ("Before HiCMA Check Rank" "\tPxQ= %3d %-3d NB= %4d N= %7d maxrank= %d "
                           "band_size= %d send_full= %d\n\n",
                           P, Q, NB, N, maxrank, band_size, send_full_tile));
#endif

#if GATHER_RANK
    /* If band_size changed, gather rank before Cholesky */
    if( auto_band && band_size > 1 ) {
        imaxrk = -1;
        iminrk=-1;
        iavgrk = -1.0;
        tlr_rank_stat("init_rank_tile", rank_array, band_size, parsec, &dcAr, &iminrk, &imaxrk, &iavgrk, maxrank);
    }
#endif

    double time_hicma = 0.0;
    info = 0;
    if( check ) {
	    /* Timer start */
	    SYNC_TIME_START();

	    /* Uncompresses approximate matrix into dcA0 */
	    STARSH_check(parsec, uplo, (parsec_tiled_matrix_dc_t*)&dcA0,
			    (parsec_tiled_matrix_dc_t*)&dcA,
			    (parsec_tiled_matrix_dc_t*)&dcAr,
			    band_size, maxrank, &info);

	    /* Timer end */
	    SYNC_TIME_PRINT(rank, ("STARSH_check" "\tPxQ= %3d %-3d NB= %4d N= %7d\n\n",
				    P, Q, NB, N));

	    /* Timer start */
	    SYNC_TIME_START();

	    /* Dense Cholesky on Approximate A0 */
	    dplasma_dpotrf(parsec, uplo, (parsec_tiled_matrix_dc_t*)&dcA0); 

	    /* Timer end */
	    SYNC_TIME_PRINT(rank, ("Dense_potrf" "\tPxQ= %3d %-3d NB= %4d N= %7d : %14f gflops\n",
				    P, Q, NB, N, gflops=(flops/1e9)/sync_time_elapsed));
    }

    if( 0 == rank && verbose ) { printf("%d: HiCMA dpotrf is starting\n", __LINE__);fflush(stdout);}

    /* Gather time info during Cholesky */
    double *g_time = calloc(6, sizeof(double));
    double *critical_path_time = calloc(6, sizeof(double));

#ifdef PARSEC_STATS_SCHED
    parsec_sched_stat_reset(parsec);
#endif
#ifdef PARSEC_STATS_COMM
    parsec_comm_stat_reset(parsec);
#endif

    /* Timer start */
    SYNC_TIME_START();

    /* HiCMA Cholesky */
    info = HiCMA_dpotrf_L(parsec, uplo, (parsec_tiled_matrix_dc_t*)&dcA,
		    (parsec_tiled_matrix_dc_t*)&dcAr,
		    (parsec_tiled_matrix_dc_t*)&dcRank,
		    tol, fixedrk, maxrank, &lookahead, band_size,
		    HNB, compmaxrank, send_full_tile, &two_flow,
		    tileopcounters, opcounters, critical_path_time, verbose
		    );

    /* Timer end */
    SYNC_TIME_PRINT(rank, ("HiCMA_potrf" "\tsend_full_tile= %d band_size= %d lookahead= %d "
			    "kind_of_problem= %d HNB= %d PxQ= %3d %-3d "
			    "NB= %4d N= %7d 2flow= %d : %14f gflops\n",
			    send_full_tile, band_size, lookahead, kind_of_problem,
			    HNB, P, Q, NB, N, two_flow, gflops=(flops/1e9)/sync_time_elapsed));
    /* Record time */
    time_hicma = sync_time_elapsed;

#ifdef PARSEC_STATS_SCHED
    parsec_sched_stat_print(parsec);
#endif
#ifdef PARSEC_STATS_COMM
    parsec_comm_stat_print(parsec);
#endif

    if( 0 == rank && info != 0 ) {
            printf("-- Factorization is suspicious (info = %d) ! \n", info);
            ret |= 1;
    }

#if PRINT_RANK
    /* Timer start */
    SYNC_TIME_START();

    parsec_rank_print(parsec, (parsec_tiled_matrix_dc_t*)&dcRank, band_size);

    /* Timer end */
    SYNC_TIME_PRINT(rank, ("rank_print" "\tPxQ= %3d %-3d NB= %4d N= %7d maxrank= %d "
			    "send_full= %d band_size= %d\n\n",
			    P, Q, NB, N, maxrank, send_full_tile, band_size));
#endif

#if DEBUG_INFO
    /* Timer start */
    SYNC_TIME_START();

    /* Check correctness of rank */
    parsec_rank_check(parsec, (parsec_tiled_matrix_dc_t*)&dcAr, band_size); 

    /* Timer end */
    SYNC_TIME_PRINT(rank, ("After HiCMA Check Rank" "\tPxQ= %3d %-3d NB= %4d N= %7d maxrank= %d "
			    "band_size= %d send_full= %d\n\n",
			    P, Q, NB, N, maxrank, band_size, send_full_tile));
#endif

    if( check ) {
	    if(rank == 0) printf("\nCHECK RESULT:\n\n");

	    /* Equality of two matrices coming from dense and low rank potrf
	    */
	    /* Uncompresses approximate matrix into dcAd */
	    if(rank == 0) printf("Uncompress dcAd\n");

	    STARSH_check(parsec, uplo, (parsec_tiled_matrix_dc_t*)&dcAd,
			    (parsec_tiled_matrix_dc_t*)&dcA,
			    (parsec_tiled_matrix_dc_t*)&dcAr,
			    band_size, maxrank, &info);

	    if(rank == 0) {
		    printf("+++++++++++++++++++++++++++\n");
            printf("Checking HiCMA L vs dense L\n");
            printf("dcAd is uncompressed matrix L obtained from HiCMA_dpotrf_L.jdf on approximate A.\n");
            printf("dcA0  is L factor obtained from dense cholesky (dpotrf_L.jdf in this repo) on approximate A.\n");
            printf("difference between these two matrices must be less than fixed accuracy threshold %.1e if fixed rank %d is zero.\n", tol, fixedrk);
        }

        ret |= check_diff( parsec, (rank == 0) ? verbose : 0, uplo,
                             (parsec_tiled_matrix_dc_t *)&dcAd,
                             (parsec_tiled_matrix_dc_t *)&dcA0, tol);
    }

    if( !info && check ) { 
        /* Check the factorization 
         * Input: initial A and computed L where A=LLT
         * Performs A-L*L^T
         * Reports norms
         */
        /* Check the factorization obtained from HiCMA */
        sym_two_dim_block_cyclic_t dcA2;
        sym_two_dim_block_cyclic_init(&dcA2, matrix_RealDouble,
                                      nodes, rank, NB, NB, N, N, 0, 0,
                                      N, N, P, uplo);
        dcA2.mat = parsec_data_allocate((size_t)dcA2.super.nb_local_tiles *
                                        (size_t)dcA2.super.bsiz *
                                        (size_t)parsec_datadist_getsizeoftype(dcA2.super.mtype));
        parsec_data_collection_set_key((parsec_data_collection_t*)&dcA2, "dcA2");

        /* Generates original problem withOUT approximation in dcA2 */
        starsh_generate_map(parsec, uplo, (parsec_tiled_matrix_dc_t *)&dcA2,
                &params);

        if(rank == 0) {
            printf("++++++++++++++++++++++++++\n");
            printf("Checking HiCMA L vs original A via computing LLT\n");
            printf("dcA2 is original problem A withOUT any approximation.\n");
            printf("dcAd is uncompressed matrix L obtained from HiCMA_dpotrf_L.jdf on approximate A.\n");
            printf("difference between these two matrices must be less than fixed accuracy threshold %.1e if fixed rank %d is zero.\n", tol, fixedrk);
        }

        /** dcAd must be filled again because it is changed somehow above */
        STARSH_check(parsec, uplo, (parsec_tiled_matrix_dc_t*)&dcAd,
                             (parsec_tiled_matrix_dc_t*)&dcA,
                             (parsec_tiled_matrix_dc_t*)&dcAr,
                             band_size, maxrank, &info);

        /* ~/parsec/build/dplasma/lib/dplasma_dcheck.c */
        ret |= check_dpotrf2( parsec, (rank == 0) ? verbose : 0, uplo,
                             (parsec_tiled_matrix_dc_t *)&dcAd, //hicma_potrf(A approximate)
                             (parsec_tiled_matrix_dc_t *)&dcA2, tol);

        parsec_data_free(dcA0.mat);
        parsec_tiled_matrix_dc_destroy( (parsec_tiled_matrix_dc_t*)&dcA0);

        parsec_data_free(dcA2.mat);
        parsec_tiled_matrix_dc_destroy( (parsec_tiled_matrix_dc_t*)&dcA2);

        parsec_data_free(dcAd.mat);
        parsec_tiled_matrix_dc_destroy( (parsec_tiled_matrix_dc_t*)&dcAd);
    }

    int fmaxrk = -1, fminrk=-1;
    double favgrk = -1.0;
#if GATHER_RANK
    tlr_rank_stat("rank_tile", rank_array, band_size, parsec, &dcAr, &fminrk, &fmaxrk, &favgrk, maxrank);
#endif

    /* Count operations */
    SYNC_TIME_START();
    unsigned long* alltileopcounters = calloc(nelm_tileopcounters, sizeof(unsigned long));
#if   defined(PARSEC_HAVE_MPI)
    MPI_Reduce(tileopcounters, alltileopcounters, nelm_tileopcounters, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);  
#elif defined(PARSEC_HAVE_LCI)
    memcpy(alltileopcounters, tileopcounters, nelm_tileopcounters * sizeof(unsigned long));
    lci_allreducel(alltileopcounters, nelm_tileopcounters * sizeof(unsigned long), lci_sum_op_ul);
#else
    memcpy(alltileopcounters, tileopcounters, nelm_tileopcounters * sizeof(unsigned long));
#endif
    
    unsigned long* allopcounters     = calloc(nelm_opcounters,     sizeof(unsigned long));
#if   defined(PARSEC_HAVE_MPI)
    MPI_Reduce(opcounters, allopcounters, nelm_opcounters, MPI_UNSIGNED_LONG, MPI_MAX, 0, MPI_COMM_WORLD);  
#elif defined(PARSEC_HAVE_LCI)
    memcpy(allopcounters, opcounters, nelm_opcounters * sizeof(unsigned long));
    lci_allreducel(allopcounters, nelm_opcounters * sizeof(unsigned long), lci_max_op_ul);
#else
    memcpy(allopcounters, opcounters, nelm_opcounters * sizeof(unsigned long));
#endif

    /* Parameters for band and critcal path operation count */
    unsigned long sum_band = 0UL;
    unsigned long total_band = 0UL;
    unsigned long sum_offband = 0UL;
    unsigned long total_offband = 0UL;
    unsigned long total_path = 0UL;
    unsigned long total_offpath = 0UL;

    /* Sum ops in a precess */
    for( int i = 1; i < cores; i++) {
        op_band[0] += op_band[i]; 
        op_offband[0] += op_offband[i]; 
        op_path[0] += op_path[i]; 
        op_offpath[0] += op_offpath[i]; 
    }

    /* Reduce to process 0 */
#if   defined(PARSEC_HAVE_MPI)
    MPI_Reduce(&op_band[0], &total_band, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&op_offband[0], &total_offband, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&op_path[0], &total_path, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&op_offpath[0], &total_offpath, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
#elif defined(PARSEC_HAVE_LCI)
    total_band    = op_band[0];
    total_offband = op_offband[0];
    total_path    = op_path[0];
    total_offpath = op_offpath[0];
    lci_allreducem(&total_band,    sizeof(unsigned long), lci_sum_op_ul);
    lci_allreducem(&total_offband, sizeof(unsigned long), lci_sum_op_ul);
    lci_allreducem(&total_path,    sizeof(unsigned long), lci_sum_op_ul);
    lci_allreducem(&total_offpath, sizeof(unsigned long), lci_sum_op_ul);
#else
    total_band = op_band[0];
    total_offband = op_offband[0];
    total_path = op_path[0];
    total_offpath = op_offpath[0];
#endif

    /* Synchronization */ 
#if   defined(PARSEC_HAVE_MPI)
    MPI_Barrier(MPI_COMM_WORLD);
#elif defined(PARSEC_HAVE_LCI)
    lci_barrier();
#endif

    unsigned long total_numop = 0;
    if(rank == 0) {
        for(int i = 0; i < NT; i++){
            for(int j = 0; j < NT; j++){
                total_numop += alltileopcounters[i*NT+j];
            }
        }

        unsigned long sum_thread = 0;
        for(int i = 0; i < nodes; i++){
            for(int j = 0; j < cores; j++){
                sum_thread += allopcounters[i*cores+j];
            }
        }

        assert(total_numop == sum_thread);

        /* operation count for tiles on band */ 
        for(int i = 0; i < NT; i++) {
            for(int j = my_own_max(i-band_size+1, 0); j <= i ; j++) {
                sum_band += alltileopcounters[i*NT+j];
            }
        }

        /* operation count for tiles off band */
        for(int i = band_size; i < NT; i++) {
            for(int j = 0; j <= i-band_size ; j++) {
                sum_offband += alltileopcounters[i*NT+j];
            }
        }

        assert(sum_band + sum_offband == total_numop);
        assert(sum_band == total_band);
        assert(sum_offband == total_offband);
        assert(total_path + total_offpath == total_numop);

    }

    SYNC_TIME_PRINT(rank, ("Operation_counts : band_size= %d total= %le op_band= %le "
                           "op_offband= %le ratio_band= %lf op_critical_path= %le "
                           "op_off_critical_path= %le ratio_critical_path= %lf\n",
                           band_size, (double)total_numop, (double)total_band, (double)total_offband, (double)total_band/total_offband,
                           (double)total_path, (double)total_offpath, (double)total_path/total_offpath));

    /* Time for critical path */
    SYNC_TIME_START();
#if   defined(PARSEC_HAVE_MPI)
    MPI_Reduce(critical_path_time, g_time, 6, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#elif defined(PARSEC_HAVE_LCI)
    memcpy(g_time, critical_path_time, 6 * sizeof(double));
    lci_allreducem(g_time, 6 * sizeof(double), lci_sum_op_d);
#else
    memcpy(g_time, critical_path_time, 6 * sizeof(double));
#endif
    SYNC_TIME_PRINT(rank, ("Total_critical_path_time %lf, potrf %lf, trsm %lf, syrk %lf\n\n", g_time[0] + g_time[1] + g_time[2], g_time[0], g_time[1], g_time[2]));

    /* Print info during computation */
    if(rank == 0) {
        printf("R-LO ");
        printf("%d %d %d %d   ", N, N, NB, NB);
        printf("%d %d   ", HNB, HNB);
        printf("%d %d %d %d   ", nodes, P, Q, cores);
        printf("%d %s %.1e %g %d %g   ", kind_of_problem, str_problem[kind_of_problem], tol, add_diag, fixedrk, wave_k);
        printf("%d %d %d %d %d ", send_full_tile, band_size, lookahead, 0, auto_band);
        printf("%d %d %d   ", maxrank, genmaxrank, compmaxrank);
        printf("%g %d %d   ", iavgrk_before, iminrk_before, imaxrk_before);
        printf("%g %d %d   ", iavgrk, iminrk, imaxrk);
        printf("%g %d %d   ", favgrk, fminrk, fmaxrk);
        printf("%lf %lf   ", memory_per_node, memory_per_node_max);
        printf("%lf %lf %lf %lf    ", g_time[0] + g_time[1] + g_time[2], g_time[0], g_time[1], g_time[2]);
        printf("%lf %lf %lf %lf %lf    ", time_starsh, time_hicma, time_opt_band, time_regenerate, 0.0);
        printf("%le %le %le %le %le ", (double)total_numop, (double)total_band, (double)total_offband, (double)total_path, (double)total_offpath);
        printf("%lf %lf %d  ", (double)total_band/total_offband, (double)total_path/total_offpath, two_flow);
       if (kind_of_problem== 6 || kind_of_problem==7){
        printf("%d %d %d  ", order, numobj, rbf_kernel);
        printf("%g %g ", radius, density);
        }
#ifdef GITHASH
        printf("%s ", xstr(GITHASH));
#else
        printf("GITHASH:N/A    ");
#endif
        printf("\"");
        for(int i=0; i<argc; i++){
            printf("%s ",  argv[i]);
        }
        printf("\"");
        printf("\n");
    }

    /* Free memory */
    free(op_band);
    free(op_offband);
    free(op_path);
    free(op_offpath);
    free(alltileopcounters);
    free(allopcounters);
    free(rank_array);
    free(g_time);
    free(critical_path_time);

#if BAND_MEMORY_CONTIGUOUS
    parsec_data_free(dcA.band.mat);
#else
    parsec_band_free(parsec, (parsec_tiled_matrix_dc_t *)&dcA, band_size, indicator_band);
#endif

    parsec_data_free(dcAr.mat);
    parsec_band_free(parsec, (parsec_tiled_matrix_dc_t *)&dcA, band_size, indicator_offband);
    parsec_band_free(parsec, (parsec_tiled_matrix_dc_t *)&dcRank, band_size, indicator_offband);
    parsec_tiled_matrix_dc_destroy( (parsec_tiled_matrix_dc_t*)&dcA );
    parsec_tiled_matrix_dc_destroy( &dcA.band.super );
    parsec_tiled_matrix_dc_destroy( (parsec_tiled_matrix_dc_t*)&dcRank );
    parsec_tiled_matrix_dc_destroy( &dcRank.band.super );
    parsec_tiled_matrix_dc_destroy( (parsec_tiled_matrix_dc_t*)&dcAr );

    cleanup_parsec(parsec, iparam);
    // Free STARS-H data
    free(index);
    if(kind_of_problem == 0)
    {
        starsh_randtlr_free(data);
    }
    else if(kind_of_problem == 1 || kind_of_problem == 5)
    {
        starsh_eddata_free(data);
    }
    else if(kind_of_problem == 2 || kind_of_problem == 3 || kind_of_problem == 4)
    {
        starsh_ssdata_free(data);
    }
    else if(kind_of_problem == 6 || kind_of_problem == 7)
    {
        starsh_mddata_free(data);
    }
    return ret;
}
