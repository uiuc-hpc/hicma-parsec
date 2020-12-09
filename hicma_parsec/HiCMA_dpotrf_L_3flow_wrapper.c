/**
 * @copyright (c) 2020 King Abdullah University of Science and Technology (KAUST).
 *                     The Universiy of Tennessee and The Universiy of Tennessee Research Foundation.
 *                     All rights reserved.
 **/

#include "hicma_parsec.h"
#include "HiCMA_dpotrf_L_3flow.h"

/* Used for gathering time during execution */
extern double *gather_time;
extern double *gather_time_tmp;

/* Task clase ID */
static int potrf_id = 3;
static int trsm_id = 4;
static int syrk_id = 5;
static int gemm_id = 6;

/* Wrap function and wrap complete function
 * to calculate time for each task:
 * potrf, trsm, syrk and gemm
 */
static int wrap_potrf(parsec_execution_stream_t * es, 
                      __parsec_HiCMA_dpotrf_L_3flow_potrf_dpotrf_task_t * this_task)
{
    parsec_HiCMA_dpotrf_L_3flow_taskpool_t *parsec_tp = (parsec_HiCMA_dpotrf_L_3flow_taskpool_t*)this_task->taskpool;
    /* Record start time of potrf */
    parsec_tp->_g_potrf_time_temp[0] = MPI_Wtime();
    gather_time_tmp[es->th_id] = parsec_tp->_g_potrf_time_temp[0];
    return parsec_tp->_g_wrap_potrf(es, (parsec_task_t *)this_task);
}

static int wrap_potrf_complete(parsec_execution_stream_t * es,
                      __parsec_HiCMA_dpotrf_L_3flow_potrf_dpotrf_task_t * this_task)
{
    int val;
    double end_time;
    parsec_HiCMA_dpotrf_L_3flow_taskpool_t *parsec_tp = (parsec_HiCMA_dpotrf_L_3flow_taskpool_t*)this_task->taskpool;
    val = parsec_tp->_g_wrap_potrf_complete(es, (parsec_task_t *)this_task);
    end_time = MPI_Wtime();
    parsec_tp->_g_potrf_time[0] += end_time - parsec_tp->_g_potrf_time_temp[0]; 
    gather_time[es->th_id] += end_time - gather_time_tmp[es->th_id];

#if PRINT_CRITICAL_PATH_TIME
    fprintf(stderr, "OUT_critical_path_time band_size %d Nodes %d Matrix %d POTRF %d end_time %lf start_time %lf exe_time %lf sum_time %lf\n",
		    parsec_tp->_g_band_size, parsec_tp->_g_descAg->super.nodes, parsec_tp->_g_descAg->lm,
		    this_task->locals.k.value, end_time, parsec_tp->_g_potrf_time_temp[0],
		    end_time - parsec_tp->_g_potrf_time_temp[0], parsec_tp->_g_potrf_time[0]);
#endif

    return val;
}

static int wrap_trsm(parsec_execution_stream_t * es,
                      __parsec_HiCMA_dpotrf_L_3flow_potrf_dtrsm_task_t * this_task)
{
    parsec_HiCMA_dpotrf_L_3flow_taskpool_t *parsec_tp = (parsec_HiCMA_dpotrf_L_3flow_taskpool_t*)this_task->taskpool;
    gather_time_tmp[es->th_id] = MPI_Wtime();

    if(this_task->locals.m.value == this_task->locals.k.value + 1){
        /* Record start time of trsm */
        parsec_tp->_g_trsm_time_temp[0] = MPI_Wtime();
        return parsec_tp->_g_wrap_trsm(es, (parsec_task_t *)this_task);
    } else {
        return parsec_tp->_g_wrap_trsm(es, (parsec_task_t *)this_task);
    }
}

static int wrap_trsm_complete(parsec_execution_stream_t * es,
                      __parsec_HiCMA_dpotrf_L_3flow_potrf_dtrsm_task_t * this_task)
{
    parsec_HiCMA_dpotrf_L_3flow_taskpool_t *parsec_tp = (parsec_HiCMA_dpotrf_L_3flow_taskpool_t*)this_task->taskpool;
    int val;
    double end_time;

    if(this_task->locals.m.value == this_task->locals.k.value + 1){
        val = parsec_tp->_g_wrap_trsm_complete(es, (parsec_task_t *)this_task);
        end_time = MPI_Wtime();
        parsec_tp->_g_trsm_time[0] += end_time - parsec_tp->_g_trsm_time_temp[0];

#if PRINT_CRITICAL_PATH_TIME
	fprintf(stderr, "OUT_critical_path_time band_size %d Nodes %d Matrix %d TRSM %d end_time %lf start_time %lf exe_time %lf sum_time %lf\n",
			parsec_tp->_g_band_size, parsec_tp->_g_descAg->super.nodes, parsec_tp->_g_descAg->lm,
			this_task->locals.k.value, end_time, parsec_tp->_g_trsm_time_temp[0],
			end_time - parsec_tp->_g_trsm_time_temp[0], parsec_tp->_g_trsm_time[0]);
#endif
    } else {
        val = parsec_tp->_g_wrap_trsm_complete(es, (parsec_task_t *)this_task);
        end_time = MPI_Wtime();
    }

    gather_time[es->th_id] += end_time - gather_time_tmp[es->th_id];
    return val;
}

static int wrap_syrk(parsec_execution_stream_t * es,
                      __parsec_HiCMA_dpotrf_L_3flow_potrf_dsyrk_task_t * this_task)
{
    parsec_HiCMA_dpotrf_L_3flow_taskpool_t *parsec_tp = (parsec_HiCMA_dpotrf_L_3flow_taskpool_t*)this_task->taskpool;
    gather_time_tmp[es->th_id] = MPI_Wtime();
    if(this_task->locals.m.value == this_task->locals.k.value + 1){
        /* Record start time of syrk */
        parsec_tp->_g_syrk_time_temp[0] = MPI_Wtime();
        return parsec_tp->_g_wrap_syrk(es, (parsec_task_t *)this_task);
    } else {
        return parsec_tp->_g_wrap_syrk(es, (parsec_task_t *)this_task);
    }
}

static int wrap_syrk_complete(parsec_execution_stream_t * es,
                      __parsec_HiCMA_dpotrf_L_3flow_potrf_dsyrk_task_t * this_task)
{
    parsec_HiCMA_dpotrf_L_3flow_taskpool_t *parsec_tp = (parsec_HiCMA_dpotrf_L_3flow_taskpool_t*)this_task->taskpool;
    int val;
    double end_time;

    if(this_task->locals.m.value == this_task->locals.k.value + 1){
        val = parsec_tp->_g_wrap_syrk_complete(es, (parsec_task_t *)this_task);
        end_time = MPI_Wtime();
        parsec_tp->_g_syrk_time[0] += end_time - parsec_tp->_g_syrk_time_temp[0];

#if PRINT_CRITICAL_PATH_TIME
	fprintf(stderr, "OUT_critical_path_time band_size %d Nodes %d Matrix %d SYRK %d end_time %lf start_time %lf exe_time %lf sum_time %lf\n",
			parsec_tp->_g_band_size, parsec_tp->_g_descAg->super.nodes, parsec_tp->_g_descAg->lm,
			this_task->locals.k.value, end_time, parsec_tp->_g_syrk_time_temp[0],
			end_time - parsec_tp->_g_syrk_time_temp[0], parsec_tp->_g_syrk_time[0]);
#endif
    } else {
        val = parsec_tp->_g_wrap_syrk_complete(es, (parsec_task_t *)this_task);
        end_time = MPI_Wtime();
    }

    gather_time[es->th_id] += end_time - gather_time_tmp[es->th_id];
    return val;
}

static int wrap_gemm(parsec_execution_stream_t * es,
                      __parsec_HiCMA_dpotrf_L_3flow_potrf_dgemm_task_t * this_task)
{
    parsec_HiCMA_dpotrf_L_3flow_taskpool_t *parsec_tp = (parsec_HiCMA_dpotrf_L_3flow_taskpool_t*)this_task->taskpool;
    /* Record start time of gemm */
    gather_time_tmp[es->th_id] = MPI_Wtime();
    return parsec_tp->_g_wrap_gemm(es, (parsec_task_t *)this_task);
}

static int wrap_gemm_complete(parsec_execution_stream_t * es,
                      __parsec_HiCMA_dpotrf_L_3flow_potrf_dgemm_task_t * this_task)
{
    parsec_HiCMA_dpotrf_L_3flow_taskpool_t *parsec_tp = (parsec_HiCMA_dpotrf_L_3flow_taskpool_t*)this_task->taskpool;
    int val;
    double start_time = gather_time_tmp[es->th_id];
    val = parsec_tp->_g_wrap_gemm_complete(es, (parsec_task_t *)this_task);
    double end_time = MPI_Wtime();
    gather_time[es->th_id] += end_time - start_time; 

    if( DEBUG_INFO )
        fprintf(stderr, "band_size %d Nodes %d Matrix %d GEMM %d %d %d end_time %lf start_time %lf exe_time %lf sum_time %lf\n",
			parsec_tp->_g_band_size, parsec_tp->_g_descAg->super.nodes, parsec_tp->_g_descAg->lm,
                        this_task->locals.m.value, this_task->locals.n.value, this_task->locals.k.value,
                        end_time, start_time, end_time - start_time, gather_time[es->th_id]);
    return val;
}

#if defined(PARSEC_HAVE_CUDA)

/* Allocate memory for workspace */
static parsec_potrf_workspace_t* workspace_memory_allocate( parsec_potrf_workspace_t *ws ) {
    ws = (parsec_potrf_workspace_t *)malloc( sizeof(parsec_potrf_workspace_t) );
    ws->gpu_workspace = (parsec_potrf_gpu_workspace_t *)malloc( parsec_nb_devices * sizeof(parsec_potrf_gpu_workspace_t) );
    
    for( int i = 0; i < parsec_nb_devices; i++ ) {
        ws->gpu_workspace[i].stream_workspace = (parsec_potrf_stream_workspace_t *)malloc( PARSEC_MAX_STREAMS * sizeof(parsec_potrf_stream_workspace_t) );
        ws->gpu_workspace[i].gpu_device = (parsec_device_cuda_module_t *)malloc( sizeof(parsec_device_cuda_module_t) ); 
    }

    return ws;
}

/* Shared with cusolver_potrf_ws */
static void workspace_memory_free( parsec_potrf_workspace_t *ws)
{
    for(int i = 0; i < (int)parsec_nb_devices; i++) {
        parsec_device_module_t *device = parsec_mca_device_get(i);
        if( NULL == device ) continue;
        if( device->type != PARSEC_DEV_CUDA ) continue;

        for(int j = 0; j < ws->gpu_workspace[i].gpu_device->max_exec_streams; j++) {
            /* j 0, h2d; j 1, d2h */
            if( j <= 1 ) continue;

            /* Free GPU buffer */
            if( NULL != ws->gpu_workspace[i].stream_workspace[j].gpu_buffer ) {
                zone_free( ws->gpu_workspace[i].gpu_device->memory, ws->gpu_workspace[i].stream_workspace[j].gpu_buffer );
            }

            /* Free GPU handle */
            if( NULL != ws->gpu_workspace[i].stream_workspace[j].handle ) {
                cusolverDnHandle_t handle = ws->gpu_workspace[i].stream_workspace[j].handle;
                cusolverStatus_t status = cusolverDnDestroy(handle);
                assert(status == CUSOLVER_STATUS_SUCCESS);
            }
        }
    }

    for( int i = 0; i < parsec_nb_devices; i++ ) {
        free( ws->gpu_workspace[i].stream_workspace ); 
    }

    free( ws->gpu_workspace );
    free( ws );
}
#endif /* PARSEC_HAVE_CUDA */

/* 3flow version, non-blocking 
 *
 * @param [in] parsec:             parsec context 
 * @param [inout] A:               the data, already distributed and allocated
 * @param [inout] Ar:              the rank info, already distributed and allocated
 * @param [in] RG:                 used for reordering GEMM 
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
HiCMA_dpotrf_L_3flow_New( parsec_context_t *parsec,
                  int uplo,
                  parsec_tiled_matrix_dc_t *A,
                  parsec_tiled_matrix_dc_t *Au,
                  parsec_tiled_matrix_dc_t *Av,
                  parsec_tiled_matrix_dc_t *Ar,
                  parsec_tiled_matrix_dc_t *RG,
                  parsec_tiled_matrix_dc_t *Rank,
                  double acc, /* accuracy threshold */
                  int fixed_rank,     /* rank threshold     */
                  int storagemaxrank,  /* size of storage    */
                  int lookahead,
                  int band_size,
                  int hmb,
                  int compmaxrank, /* size of temporary buffers */
                  int send_full_tile,
                  unsigned long* tileopcounters,
                  unsigned long* opcounters,
                  double *critical_path_time,
                  int *info )
{
    parsec_taskpool_t *tp = NULL;
    void** hook;

    /* Check input arguments */
    if ((uplo != PlasmaUpper) && (uplo != PlasmaLower)) {
        dplasma_error("HiCMA_dpotrf_L_3flow_New", "illegal value of uplo");
        return NULL /*-1*/;
    }
    if ((uplo == PlasmaUpper)) {
        dplasma_error("HiCMA_dpotrf_L_3flow_New", "HiCMA_dpotrf_L_3flow does not support PlasmaUpper for now");
        return NULL /*-1*/;
    }

    if(storagemaxrank > compmaxrank) {
        dplasma_error("HiCMA_dpotrf_L_3flow_New", "maxrank for storage larger than maxrank for buffers \
                is not meaningful");
        return NULL /*-1*/;
    }
    if(storagemaxrank > (A->mb/2)) {
        dplasma_error("HiCMA_dpotrf_L_3flow_New", "maxrank should not be larger than half of block size");
        return NULL /*-1*/;
    }

    int nb = 0;
#if defined(PARSEC_HAVE_CUDA)
    /** Find all CUDA devices */
    for(int i = 0; i < (int)parsec_nb_devices; i++) {
        parsec_device_module_t *device = parsec_mca_device_get(i);
        if( PARSEC_DEV_CUDA == device->type ) {
            nb++;
        }
    }
    if(nb == 0) {
        char hostname[256];
        gethostname(hostname, 256);
        fprintf(stderr, "No CUDA device found on rank %d on %s\n",
                parsec->my_rank, hostname);
    }
    int *dev_index = (int*)malloc(nb * sizeof(int));
    nb = 0;
    for(int i = 0; i < (int)parsec_nb_devices; i++) {
        parsec_device_module_t *device = parsec_mca_device_get(i);
        if( PARSEC_DEV_CUDA == device->type ) {
            dev_index[nb++] = device->device_index;
        }
    }

    /* Declare workspace used on GPU */
    parsec_potrf_workspace_t *ws_handle, *ws_mbr, *ws_rr;

    /* Allocate memory */
    ws_handle = workspace_memory_allocate( ws_handle ); 
    ws_mbr = workspace_memory_allocate( ws_mbr ); 
    ws_rr = workspace_memory_allocate( ws_rr ); 

    /* Traverse all gpu device */
    for(int i = 0; i < (int)parsec_nb_devices; i++) {
        parsec_device_module_t *device = parsec_mca_device_get(i);
        if( NULL == device ) continue;
        if( device->type != PARSEC_DEV_CUDA ) continue;

        parsec_device_cuda_module_t *gpu_device = (parsec_device_cuda_module_t*)device;
        cudaSetDevice(gpu_device->cuda_index);

        ws_handle->gpu_workspace[i].gpu_device = gpu_device;
        ws_mbr->gpu_workspace[i].gpu_device = gpu_device;
        ws_rr->gpu_workspace[i].gpu_device = gpu_device;

        /* Traverse all streams */ 
        for(int j = 0; j < gpu_device->max_exec_streams; j++) {
            /* j 0, h2d; j 1, d2h */
            if( j <= 1 ) continue;

            cublasStatus_t status;
            cudaError_t cudaStatus;
            cusolverDnHandle_t handle;

            /* Create handle */
            {
                status = cusolverDnCreate(&handle);
		assert(CUSOLVER_STATUS_SUCCESS == status);
		status = cusolverDnSetStream(handle, gpu_device->exec_stream[j].cuda_stream);
		assert(CUSOLVER_STATUS_SUCCESS == status);
		ws_handle->gpu_workspace[i].stream_workspace[j].handle = handle;
		ws_mbr->gpu_workspace[i].stream_workspace[j].handle = NULL;
		ws_rr->gpu_workspace[i].stream_workspace[j].handle = NULL;
            }

            /* Allocate workspace for potrf handle */
            {
                int workspace_size;
                status = cusolverDnDpotrf_bufferSize(handle, CUBLAS_FILL_MODE_LOWER, A->nb, NULL, A->mb, &workspace_size);
                assert(CUSOLVER_STATUS_SUCCESS == status);
                ws_handle->gpu_workspace[i].stream_workspace[j].gpu_buffer = zone_malloc( gpu_device->memory, workspace_size * sizeof(double) + sizeof(int) );
                assert(NULL != ws_handle->gpu_workspace[i].stream_workspace[j].gpu_buffer);
                ws_handle->gpu_workspace[i].stream_workspace[j].buffer_size = workspace_size;
            }

            /* Temporary buffer for size mb * compmaxrank */
            {
                int workspace_size = A->mb * compmaxrank;
                ws_mbr->gpu_workspace[i].stream_workspace[j].gpu_buffer = zone_malloc( gpu_device->memory, workspace_size * sizeof(double) + sizeof(int) );
                assert(NULL != ws_mbr->gpu_workspace[i].stream_workspace[j].gpu_buffer);
                ws_mbr->gpu_workspace[i].stream_workspace[j].buffer_size = workspace_size; 
            }

            /* Temporary buffer for size compmaxrank * compmaxrank */
            {
                int workspace_size = compmaxrank * compmaxrank;
                ws_rr->gpu_workspace[i].stream_workspace[j].gpu_buffer = zone_malloc( gpu_device->memory, workspace_size * sizeof(double) + sizeof(int) );
                assert(NULL != ws_rr->gpu_workspace[i].stream_workspace[j].gpu_buffer);
                ws_rr->gpu_workspace[i].stream_workspace[j].buffer_size = workspace_size; 
            }
        }
    }
#endif

    *info = 0;
    if ( uplo == PlasmaLower ) {
        parsec_HiCMA_dpotrf_L_3flow_taskpool_t *hicma_dpotrf =
            parsec_HiCMA_dpotrf_L_3flow_new( uplo, A, Au, Av, Ar, Rank, info,
                                       acc, fixed_rank, storagemaxrank,
				       compmaxrank, lookahead, band_size);

#if defined(PARSEC_HAVE_CUDA)
        hicma_dpotrf->_g_ws_handle = ws_handle;
        hicma_dpotrf->_g_ws_mbr = ws_mbr;
        hicma_dpotrf->_g_ws_rr = ws_rr;
        hicma_dpotrf->_g_nb_cuda_devices = nb;
        hicma_dpotrf->_g_cuda_device_index = dev_index;
#endif

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

        /* On GPU */
        if( nb > 0 ) {
//#if defined(PARSEC_HAVE_CUDA)
#if 0
            /* GPU wrap kernel */    
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
#endif

          /* Recursive */
        } else if( hmb < A->mb ) { 
//#if defined(PARSEC_HAVE_CUDA)
#if 0 
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
#else
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
#endif

          /* Others */
        } else {
//#if defined(PARSEC_HAVE_CUDA)
#if 0 
            hicma_dpotrf->_g_wrap_potrf = hicma_dpotrf->super.task_classes_array[potrf_id]->incarnations[2].hook;
            hicma_dpotrf->_g_wrap_trsm  = hicma_dpotrf->super.task_classes_array[trsm_id]->incarnations[2].hook;
            hicma_dpotrf->_g_wrap_syrk  = hicma_dpotrf->super.task_classes_array[syrk_id]->incarnations[2].hook;
            hicma_dpotrf->_g_wrap_gemm  = hicma_dpotrf->super.task_classes_array[gemm_id]->incarnations[2].hook;
            hook  = (void *)&hicma_dpotrf->super.task_classes_array[potrf_id]->incarnations[2].hook;
            *hook = &wrap_potrf;
            hook  = (void *)&hicma_dpotrf->super.task_classes_array[trsm_id]->incarnations[2].hook;
            *hook = &wrap_trsm;
            hook  = (void *)&hicma_dpotrf->super.task_classes_array[syrk_id]->incarnations[2].hook;
            *hook = &wrap_syrk;
            hook  = (void *)&hicma_dpotrf->super.task_classes_array[gemm_id]->incarnations[2].hook;
            *hook = &wrap_gemm;
#else
            hicma_dpotrf->_g_wrap_potrf = hicma_dpotrf->super.task_classes_array[potrf_id]->incarnations[1].hook;
            hicma_dpotrf->_g_wrap_trsm  = hicma_dpotrf->super.task_classes_array[trsm_id]->incarnations[0].hook;
            hicma_dpotrf->_g_wrap_syrk  = hicma_dpotrf->super.task_classes_array[syrk_id]->incarnations[0].hook;
            hicma_dpotrf->_g_wrap_gemm  = hicma_dpotrf->super.task_classes_array[gemm_id]->incarnations[0].hook;
            hook  = (void *)&hicma_dpotrf->super.task_classes_array[potrf_id]->incarnations[1].hook;
            *hook = &wrap_potrf;
            hook  = (void *)&hicma_dpotrf->super.task_classes_array[trsm_id]->incarnations[0].hook;
            *hook = &wrap_trsm;
            hook  = (void *)&hicma_dpotrf->super.task_classes_array[syrk_id]->incarnations[0].hook;
            *hook = &wrap_syrk;
            hook  = (void *)&hicma_dpotrf->super.task_classes_array[gemm_id]->incarnations[0].hook;
            *hook = &wrap_gemm;
#endif
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
        //parsec_private_memory_init( hicma_dpotrf->_g_p_work, A->nb * A->nb * sizeof(double) * 8 );
        /* size used for temporary buffers obtained from hicma/compute/pzpotrf.c line 96 */
        // TODO @kadir why is this?
        size_t ws_worker = 0;
        ws_worker =  //FIXME tentative size. Find exact size. I think syrk uses less memory
        //This workspace need to be fixed, not all tasks below need it nor need that much
        2 * A->mb * 2 * compmaxrank   // for copying CU and CV into temporary buffer instead of using CUV itself. There is 2*maxrk because these buffers will be used to put two U's side by side
        + 2 * A->mb       // qrtauA qrtauB
        + compmaxrank * compmaxrank    // qrb_aubut  AcolBcolT
        + 2 * A->mb * 2 * compmaxrank // newU newV
        + (2*compmaxrank) * (2*compmaxrank)    // svd_rA     _rA
        //+ maxrk * maxrk    // svd_rB     _rB   I assume that use_trmm=1 so I commented out
        //+ maxrk * maxrk    // svd_T      _T    I assume that use_trmm=1 so I commented out
        + (2*compmaxrank)          // sigma
        + (2*compmaxrank); // superb
        ;
        ws_worker *= sizeof(double); 
        parsec_private_memory_init( hicma_dpotrf->_g_p_work, ws_worker ); 

        hicma_dpotrf->_g_p_work_rr = (parsec_memory_pool_t*)malloc(sizeof(parsec_memory_pool_t));
        parsec_private_memory_init( hicma_dpotrf->_g_p_work_rr, compmaxrank * compmaxrank * sizeof(double) ); 

        hicma_dpotrf->_g_p_work_mbr = (parsec_memory_pool_t*)malloc(sizeof(parsec_memory_pool_t));
        parsec_private_memory_init( hicma_dpotrf->_g_p_work_mbr, A->mb * compmaxrank * sizeof(double) ); 

        /* Arena */
        parsec_matrix_add2arena(&hicma_dpotrf->arenas_datatypes[PARSEC_HiCMA_dpotrf_L_3flow_FULL_ARENA],
                                parsec_datatype_double_t, matrix_UpperLower,
                                1, A->mb, A->mb, A->mb,
                                PARSEC_ARENA_ALIGNMENT_SSE, -1 );

        parsec_matrix_add2arena(&hicma_dpotrf->arenas_datatypes[PARSEC_HiCMA_dpotrf_L_3flow_DEFAULT_ARENA],
                                parsec_datatype_double_t, matrix_UpperLower,
                                1, 1, 1, 1,
                                PARSEC_ARENA_ALIGNMENT_SSE, -1 );

        parsec_matrix_add2arena(&hicma_dpotrf->arenas_datatypes[PARSEC_HiCMA_dpotrf_L_3flow_UV_ARENA],
                                parsec_datatype_double_t, matrix_UpperLower,
                                1, A->mb, storagemaxrank, A->mb,
                                PARSEC_ARENA_ALIGNMENT_SSE, -1 );

        parsec_matrix_add2arena(&hicma_dpotrf->arenas_datatypes[PARSEC_HiCMA_dpotrf_L_3flow_AR_ARENA],
                                parsec_datatype_int_t, matrix_UpperLower,
                                1, 1, 1, 1,
                                PARSEC_ARENA_ALIGNMENT_SSE, -1 );
    }

    return tp;
}

/* Destructor */
void HiCMA_dpotrf_L_3flow_Destruct(parsec_taskpool_t* _tp)
{
    parsec_HiCMA_dpotrf_L_3flow_taskpool_t *tp = (parsec_HiCMA_dpotrf_L_3flow_taskpool_t*)_tp;

#if defined(PARSEC_HAVE_CUDA)
    if( tp->_g_nb_cuda_devices > 0 ) {
        workspace_memory_free( tp->_g_ws_handle );
        workspace_memory_free( tp->_g_ws_mbr );
        workspace_memory_free( tp->_g_ws_rr );

        if( NULL != tp->_g_cuda_device_index )
            free(tp->_g_cuda_device_index);
    }
#endif

    parsec_matrix_del2arena( &tp->arenas_datatypes[PARSEC_HiCMA_dpotrf_L_3flow_DEFAULT_ARENA] );
    parsec_matrix_del2arena( &tp->arenas_datatypes[PARSEC_HiCMA_dpotrf_L_3flow_FULL_ARENA] );
    parsec_matrix_del2arena( &tp->arenas_datatypes[PARSEC_HiCMA_dpotrf_L_3flow_UV_ARENA] );
    parsec_matrix_del2arena( &tp->arenas_datatypes[PARSEC_HiCMA_dpotrf_L_3flow_AR_ARENA] );
    parsec_private_memory_fini( tp->_g_p_work );
    parsec_private_memory_fini( tp->_g_p_work_mbr );
    parsec_private_memory_fini( tp->_g_p_work_rr );

    parsec_taskpool_free(_tp);
}
