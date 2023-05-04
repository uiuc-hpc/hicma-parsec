/**
 * @copyright (c) 2021     King Abdullah University of Science and Technology (KAUST).
 * @copyright (c) 2021     The Universiy of Tennessee and The Universiy of Tennessee Research Foundation.
 *                         All rights reserved.
 * @version 0.1.0
 * @date 2021-01-24
 *
 **/

#include "hicma_parsec.h"
#include "HiCMA_dpotrf_L_2flow.h"

/* Used for gathering time during execution */
extern double *gather_time;
extern double *gather_time_tmp;

/* Task clase ID */
static int potrf_id = 2;
static int trsm_id = 3;
static int syrk_id = 4;
static int gemm_id = 5;

/* Wrap function and wrap complete function
 * to calculate time for each task:
 * potrf, trsm, syrk and gemm
 */
static int wrap_potrf(parsec_execution_stream_t * es, 
                      __parsec_HiCMA_dpotrf_L_2flow_potrf_dpotrf_task_t * this_task)
{
    parsec_HiCMA_dpotrf_L_2flow_taskpool_t *parsec_tp = (parsec_HiCMA_dpotrf_L_2flow_taskpool_t*)this_task->taskpool;
    /* Record start time of potrf */
    parsec_tp->_g_potrf_time_temp[0] = Wtime();
    gather_time_tmp[es->th_id] = parsec_tp->_g_potrf_time_temp[0];
    return parsec_tp->_g_wrap_potrf(es, (parsec_task_t *)this_task);
}

static int wrap_potrf_complete(parsec_execution_stream_t * es,
                      __parsec_HiCMA_dpotrf_L_2flow_potrf_dpotrf_task_t * this_task)
{
    int val;
    double end_time;
    parsec_HiCMA_dpotrf_L_2flow_taskpool_t *parsec_tp = (parsec_HiCMA_dpotrf_L_2flow_taskpool_t*)this_task->taskpool;
    val = parsec_tp->_g_wrap_potrf_complete(es, (parsec_task_t *)this_task);
    end_time = Wtime();
    parsec_tp->_g_potrf_time[0] += end_time - parsec_tp->_g_potrf_time_temp[0]; 
    gather_time[es->th_id] += end_time - gather_time_tmp[es->th_id];

#if PRINT_CRITICAL_PATH_TIME
    fprintf(stderr, "OUT_critical_path_time band_size %d Nodes %d Matrix %d POTRF %d end_time %lf start_time %lf exe_time %lf sum_time %lf\n",
		    parsec_tp->_g_band_size, parsec_tp->_g_descA->super.nodes, parsec_tp->_g_descA->lm,
		    this_task->locals.k.value, end_time, parsec_tp->_g_potrf_time_temp[0],
		    end_time - parsec_tp->_g_potrf_time_temp[0], parsec_tp->_g_potrf_time[0]);
#endif

    return val;
}

static int wrap_trsm(parsec_execution_stream_t * es,
                      __parsec_HiCMA_dpotrf_L_2flow_potrf_dtrsm_task_t * this_task)
{
    parsec_HiCMA_dpotrf_L_2flow_taskpool_t *parsec_tp = (parsec_HiCMA_dpotrf_L_2flow_taskpool_t*)this_task->taskpool;
    gather_time_tmp[es->th_id] = Wtime();

    if(this_task->locals.m.value == this_task->locals.k.value + 1){
        /* Record start time of trsm */
        parsec_tp->_g_trsm_time_temp[0] = Wtime();
        return parsec_tp->_g_wrap_trsm(es, (parsec_task_t *)this_task);
    } else {
        return parsec_tp->_g_wrap_trsm(es, (parsec_task_t *)this_task);
    }
}

static int wrap_trsm_complete(parsec_execution_stream_t * es,
                      __parsec_HiCMA_dpotrf_L_2flow_potrf_dtrsm_task_t * this_task)
{
    parsec_HiCMA_dpotrf_L_2flow_taskpool_t *parsec_tp = (parsec_HiCMA_dpotrf_L_2flow_taskpool_t*)this_task->taskpool;
    int val;
    double end_time;

    if(this_task->locals.m.value == this_task->locals.k.value + 1){
        val = parsec_tp->_g_wrap_trsm_complete(es, (parsec_task_t *)this_task);
        end_time = Wtime();
        parsec_tp->_g_trsm_time[0] += end_time - parsec_tp->_g_trsm_time_temp[0];

#if PRINT_CRITICAL_PATH_TIME
	fprintf(stderr, "OUT_critical_path_time band_size %d Nodes %d Matrix %d TRSM %d end_time %lf start_time %lf exe_time %lf sum_time %lf\n",
			parsec_tp->_g_band_size, parsec_tp->_g_descA->super.nodes, parsec_tp->_g_descA->lm,
			this_task->locals.k.value, end_time, parsec_tp->_g_trsm_time_temp[0],
			end_time - parsec_tp->_g_trsm_time_temp[0], parsec_tp->_g_trsm_time[0]);
#endif
    } else {
        val = parsec_tp->_g_wrap_trsm_complete(es, (parsec_task_t *)this_task);
        end_time = Wtime();
    }

    gather_time[es->th_id] += end_time - gather_time_tmp[es->th_id];
    return val;
}

static int wrap_syrk(parsec_execution_stream_t * es,
                      __parsec_HiCMA_dpotrf_L_2flow_potrf_dsyrk_task_t * this_task)
{
    parsec_HiCMA_dpotrf_L_2flow_taskpool_t *parsec_tp = (parsec_HiCMA_dpotrf_L_2flow_taskpool_t*)this_task->taskpool;
    gather_time_tmp[es->th_id] = Wtime();
    if(this_task->locals.m.value == this_task->locals.k.value + 1){
        /* Record start time of syrk */
        parsec_tp->_g_syrk_time_temp[0] = Wtime();
        return parsec_tp->_g_wrap_syrk(es, (parsec_task_t *)this_task);
    } else {
        return parsec_tp->_g_wrap_syrk(es, (parsec_task_t *)this_task);
    }
}

static int wrap_syrk_complete(parsec_execution_stream_t * es,
                      __parsec_HiCMA_dpotrf_L_2flow_potrf_dsyrk_task_t * this_task)
{
    parsec_HiCMA_dpotrf_L_2flow_taskpool_t *parsec_tp = (parsec_HiCMA_dpotrf_L_2flow_taskpool_t*)this_task->taskpool;
    int val;
    double end_time;

    if(this_task->locals.m.value == this_task->locals.k.value + 1){
        val = parsec_tp->_g_wrap_syrk_complete(es, (parsec_task_t *)this_task);
        end_time = Wtime();
        parsec_tp->_g_syrk_time[0] += end_time - parsec_tp->_g_syrk_time_temp[0];

#if PRINT_CRITICAL_PATH_TIME
	fprintf(stderr, "OUT_critical_path_time band_size %d Nodes %d Matrix %d SYRK %d end_time %lf start_time %lf exe_time %lf sum_time %lf\n",
			parsec_tp->_g_band_size, parsec_tp->_g_descA->super.nodes, parsec_tp->_g_descA->lm,
			this_task->locals.k.value, end_time, parsec_tp->_g_syrk_time_temp[0],
			end_time - parsec_tp->_g_syrk_time_temp[0], parsec_tp->_g_syrk_time[0]);
#endif
    } else {
        val = parsec_tp->_g_wrap_syrk_complete(es, (parsec_task_t *)this_task);
        end_time = Wtime();
    }

    gather_time[es->th_id] += end_time - gather_time_tmp[es->th_id];
    return val;
}

static int wrap_gemm(parsec_execution_stream_t * es,
                      __parsec_HiCMA_dpotrf_L_2flow_potrf_dgemm_task_t * this_task)
{
    parsec_HiCMA_dpotrf_L_2flow_taskpool_t *parsec_tp = (parsec_HiCMA_dpotrf_L_2flow_taskpool_t*)this_task->taskpool;
    /* Record start time of gemm */
    gather_time_tmp[es->th_id] = Wtime();
    return parsec_tp->_g_wrap_gemm(es, (parsec_task_t *)this_task);
}

static int wrap_gemm_complete(parsec_execution_stream_t * es,
                      __parsec_HiCMA_dpotrf_L_2flow_potrf_dgemm_task_t * this_task)
{
    parsec_HiCMA_dpotrf_L_2flow_taskpool_t *parsec_tp = (parsec_HiCMA_dpotrf_L_2flow_taskpool_t*)this_task->taskpool;
    int val;
    double start_time = gather_time_tmp[es->th_id];
    val = parsec_tp->_g_wrap_gemm_complete(es, (parsec_task_t *)this_task);
    double end_time = Wtime();
    gather_time[es->th_id] += end_time - start_time; 

    if( DEBUG_INFO )
        fprintf(stderr, "band_size %d Nodes %d Matrix %d GEMM %d %d %d end_time %lf start_time %lf exe_time %lf sum_time %lf\n",
			parsec_tp->_g_band_size, parsec_tp->_g_descA->super.nodes, parsec_tp->_g_descA->lm,
                        this_task->locals.m.value, this_task->locals.n.value, this_task->locals.k.value,
                        end_time, start_time, end_time - start_time, gather_time[es->th_id]);
    return val;
}


/* 2flow version, non-blocking 
 *
 * @param [in] parsec:             parsec context 
 * @param [inout] A:               the data, already distributed and allocated
 * @param [inout] Ar:              the rank info, already distributed and allocated
 * @param [in] Rank:               used for gathering time info during factorization 
 * @param [in] acc:                accuracy threshold 
 * @param [in] fixed_rank:         fixed rank threshold used in recompression stage of HCORE_GEMM 
 * @param [in] maxrank:            max rank threshold for storage
 * @param [in] lookahead:          lookahead to guide workflow 
 * @param [in] band_size:          band size to control dense tiles; band_size=1: only diagonal in dense 
 * @param [in] hmb:                hierarchical mb; the sub tile size used in the recursive kernel 
 * @param [in] compmaxrank:        max rank threshold used in computation 
 * @param [in] send_full_tile:     whether send full tile during factorization; default is false to give better performance 
 * @param [in] tileopcounters:     count the number of tiles 
 * @param [in] opcounters:         count operations during factorization 
 * @param [in] critical_path_time: measure the critical path time 
 * @param [intout] info            0 on all nodes if successful. 
 *                                 > 0 if the leading minor of order i of A is not positive
 *                                 definite, so the factorization could not be completed, and the
 *                                 solution has not been computed. Info will be equal to i on the
 *                                 node that owns the diagonal element (i,i), and 0 on all other nodes
 */
parsec_taskpool_t*
HiCMA_dpotrf_L_2flow_New( parsec_context_t *parsec,
                  int uplo,
                  parsec_tiled_matrix_dc_t *A,
                  parsec_tiled_matrix_dc_t *Ar,
                  parsec_tiled_matrix_dc_t *Rank,
                  double acc,
                  int fixed_rank,
                  int storagemaxrank,
                  int lookahead,
                  int band_size,
                  int hmb,
                  int compmaxrank,
                  int send_full_tile,
                  unsigned long* tileopcounters,
                  unsigned long* opcounters,
                  double *critical_path_time,
                  int *info )
{
    parsec_taskpool_t *tp = NULL;
    void** hook;
    void** eval_gpu_potrf;
    void** eval_gpu_trsm;
    void** eval_gpu_syrk;
    void** eval_gpu_gemm;

    /* Check input arguments */
    if ((uplo != dplasmaUpper) && (uplo != dplasmaLower)) {
        dplasma_error("HiCMA_dpotrf_L_2flow_New", "illegal value of uplo");
        return NULL /*-1*/;
    }
    if ((uplo == dplasmaUpper)) {
        dplasma_error("HiCMA_dpotrf_L_2flow_New", "HiCMA_dpotrf_L_2flow does not support dplasmaUpper for now");
        return NULL /*-1*/;
    }

    if(storagemaxrank > compmaxrank) {
        dplasma_error("HiCMA_dpotrf_L_2flow_New", "maxrank for storage larger than maxrank for buffers \
                is not meaningful");
        return NULL /*-1*/;
    }

    if(storagemaxrank > (A->mb/2)) {
        dplasma_error("HiCMA_dpotrf_L_2flow_New", "maxrank should not be larger than half of block size");
        return NULL /*-1*/;
    }

    *info = 0;
    parsec_HiCMA_dpotrf_L_2flow_taskpool_t *hicma_dpotrf =
        parsec_HiCMA_dpotrf_L_2flow_new( uplo, A, Ar, Rank, info,
                acc, fixed_rank, storagemaxrank,
                compmaxrank, lookahead, band_size);

    hicma_dpotrf->_g_enable_potrf = 1;
    hicma_dpotrf->_g_enable_trsm  = 1;
    hicma_dpotrf->_g_enable_syrk  = 1;
    hicma_dpotrf->_g_enable_gemm  = 1;
    hicma_dpotrf->_g_send_full_tile = send_full_tile;
    hicma_dpotrf->_g_tileopcounters = tileopcounters;
    hicma_dpotrf->_g_opcounters = opcounters;
    hicma_dpotrf->_g_smallnb = hmb;   /* recursive tile size in po */
    tp = (parsec_taskpool_t*)hicma_dpotrf;

    /* Get critical path time */
    hicma_dpotrf->_g_potrf_time = &critical_path_time[0];
    hicma_dpotrf->_g_trsm_time = &critical_path_time[1];
    hicma_dpotrf->_g_syrk_time = &critical_path_time[2];
    hicma_dpotrf->_g_potrf_time_temp = &critical_path_time[3];
    hicma_dpotrf->_g_trsm_time_temp = &critical_path_time[4];
    hicma_dpotrf->_g_syrk_time_temp = &critical_path_time[5];

    /* Recursive */
    if( hmb < A->mb ) { 
        hicma_dpotrf->_g_wrap_potrf = hicma_dpotrf->super.task_classes_array[potrf_id]->incarnations[0].hook;
        hicma_dpotrf->_g_wrap_trsm  = hicma_dpotrf->super.task_classes_array[trsm_id]->incarnations[0].hook;
        hicma_dpotrf->_g_wrap_syrk  = hicma_dpotrf->super.task_classes_array[syrk_id]->incarnations[0].hook;
        hicma_dpotrf->_g_wrap_gemm  = hicma_dpotrf->super.task_classes_array[gemm_id]->incarnations[0].hook;
        hook  = (void *)&hicma_dpotrf->super.task_classes_array[potrf_id]->incarnations[0].hook;
        *hook = &wrap_potrf;
        hook  = (void *)&hicma_dpotrf->super.task_classes_array[trsm_id]->incarnations[0].hook;
        *hook = &wrap_trsm;
        hook  = (void *)&hicma_dpotrf->super.task_classes_array[syrk_id]->incarnations[0].hook;
        *hook = &wrap_syrk;
        hook  = (void *)&hicma_dpotrf->super.task_classes_array[gemm_id]->incarnations[0].hook;
        *hook = &wrap_gemm;

    /* Others */
    } else {

        hicma_dpotrf->_g_wrap_potrf = hicma_dpotrf->super.task_classes_array[potrf_id]->incarnations[1].hook;
        hicma_dpotrf->_g_wrap_trsm  = hicma_dpotrf->super.task_classes_array[trsm_id]->incarnations[1].hook;
        hicma_dpotrf->_g_wrap_syrk  = hicma_dpotrf->super.task_classes_array[syrk_id]->incarnations[1].hook;
        hicma_dpotrf->_g_wrap_gemm  = hicma_dpotrf->super.task_classes_array[gemm_id]->incarnations[1].hook;
        hook  = (void *)&hicma_dpotrf->super.task_classes_array[potrf_id]->incarnations[1].hook;
        *hook = &wrap_potrf;
        hook  = (void *)&hicma_dpotrf->super.task_classes_array[trsm_id]->incarnations[1].hook;
        *hook = &wrap_trsm;
        hook  = (void *)&hicma_dpotrf->super.task_classes_array[syrk_id]->incarnations[1].hook;
        *hook = &wrap_syrk;
        hook  = (void *)&hicma_dpotrf->super.task_classes_array[gemm_id]->incarnations[1].hook;
        *hook = &wrap_gemm;
    }

    hicma_dpotrf->_g_wrap_potrf_complete = hicma_dpotrf->super.task_classes_array[potrf_id]->complete_execution;
    hook  = (void *)&hicma_dpotrf->super.task_classes_array[potrf_id]->complete_execution;
    *hook = &wrap_potrf_complete;

    hicma_dpotrf->_g_wrap_trsm_complete = hicma_dpotrf->super.task_classes_array[trsm_id]->complete_execution;
    hook  = (void *)&hicma_dpotrf->super.task_classes_array[trsm_id]->complete_execution;
    *hook = &wrap_trsm_complete;

    hicma_dpotrf->_g_wrap_syrk_complete = hicma_dpotrf->super.task_classes_array[syrk_id]->complete_execution;
    hook  = (void *)&hicma_dpotrf->super.task_classes_array[syrk_id]->complete_execution;
    *hook = &wrap_syrk_complete;

    hicma_dpotrf->_g_wrap_gemm_complete = hicma_dpotrf->super.task_classes_array[gemm_id]->complete_execution;
    hook  = (void *)&hicma_dpotrf->super.task_classes_array[gemm_id]->complete_execution;
    *hook = &wrap_gemm_complete;

    /* Memory pool */
    hicma_dpotrf->_g_p_work = (parsec_memory_pool_t*)malloc(sizeof(parsec_memory_pool_t));
    /* size used for temporary buffers obtained from hicma/compute/pzpotrf.c line 96 */
    size_t ws_worker = 0;
    ws_worker = 
        2 * A->mb * 2 * compmaxrank            // for copying CU and CV into temporary buffer instead of using CUV itself. There is 2*maxrk because these buffers will be used to put two U's side by side
        + 2 * A->mb                            // qrtauA qrtauB
        + compmaxrank * compmaxrank            // qrb_aubut  AcolBcolT
        + 2 * A->mb * 2 * compmaxrank          // newU newV
        + (2*compmaxrank) * (2*compmaxrank)    // svd_rA     _rA
        + (2*compmaxrank)                      // sigma
        + (2*compmaxrank);                     // superb
    ;
    ws_worker *= sizeof(double); 
    parsec_private_memory_init( hicma_dpotrf->_g_p_work, ws_worker ); 

    hicma_dpotrf->_g_p_work_rr = (parsec_memory_pool_t*)malloc(sizeof(parsec_memory_pool_t));
    parsec_private_memory_init( hicma_dpotrf->_g_p_work_rr, compmaxrank * compmaxrank * sizeof(double) ); 

    hicma_dpotrf->_g_p_work_mbr = (parsec_memory_pool_t*)malloc(sizeof(parsec_memory_pool_t));
    parsec_private_memory_init( hicma_dpotrf->_g_p_work_mbr, A->mb * compmaxrank * sizeof(double) ); 

    hicma_dpotrf->_g_PRI_CHANGE = dplasma_aux_get_priority_limit( "POTRF", A );
    if(0 == hicma_dpotrf->_g_PRI_CHANGE) {
        hicma_dpotrf->_g_PRI_CHANGE = A->nt;
    }

    /* Arena */
    parsec_matrix_add2arena(&hicma_dpotrf->arenas_datatypes[PARSEC_HiCMA_dpotrf_L_2flow_FULL_ARENA],
            parsec_datatype_double_t, matrix_UpperLower,
            1, A->mb, A->mb, A->mb,
            PARSEC_ARENA_ALIGNMENT_SSE, -1 );

    parsec_matrix_add2arena(&hicma_dpotrf->arenas_datatypes[PARSEC_HiCMA_dpotrf_L_2flow_DEFAULT_ARENA],
            parsec_datatype_double_t, matrix_UpperLower,
            1, 1, 1, 1,
            PARSEC_ARENA_ALIGNMENT_SSE, -1 );

    parsec_matrix_add2arena(&hicma_dpotrf->arenas_datatypes[PARSEC_HiCMA_dpotrf_L_2flow_UV_ARENA],
            parsec_datatype_double_t, matrix_UpperLower,
            1, A->mb, storagemaxrank*2, A->mb,
            PARSEC_ARENA_ALIGNMENT_SSE, -1 );

    parsec_matrix_add2arena(&hicma_dpotrf->arenas_datatypes[PARSEC_HiCMA_dpotrf_L_2flow_AR_ARENA],
            parsec_datatype_int_t, matrix_UpperLower,
            1, 1, 1, 1,
            PARSEC_ARENA_ALIGNMENT_SSE, -1 );

    return tp;
}

/* Destructor */
void HiCMA_dpotrf_L_2flow_Destruct(parsec_taskpool_t* _tp)
{
    parsec_HiCMA_dpotrf_L_2flow_taskpool_t *tp = (parsec_HiCMA_dpotrf_L_2flow_taskpool_t*)_tp;

    parsec_matrix_del2arena( &tp->arenas_datatypes[PARSEC_HiCMA_dpotrf_L_2flow_DEFAULT_ARENA] );
    parsec_matrix_del2arena( &tp->arenas_datatypes[PARSEC_HiCMA_dpotrf_L_2flow_FULL_ARENA] );
    parsec_matrix_del2arena( &tp->arenas_datatypes[PARSEC_HiCMA_dpotrf_L_2flow_UV_ARENA] );
    parsec_matrix_del2arena( &tp->arenas_datatypes[PARSEC_HiCMA_dpotrf_L_2flow_AR_ARENA] );
    parsec_private_memory_fini( tp->_g_p_work );
    parsec_private_memory_fini( tp->_g_p_work_mbr );
    parsec_private_memory_fini( tp->_g_p_work_rr );

    parsec_taskpool_free(_tp);
}
