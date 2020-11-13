/*
 * Copyright (c) 2010-2018 The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 * Copyright (c) 2013      Inria. All rights reserved.
 * $COPYRIGHT
 *
 * @generated d Mon Aug 20 11:38:25 2018
 *
 */
#include "parsec/data_dist/matrix/sym_two_dim_rectangle_cyclic_band.h"
#include "parsec/data_dist/matrix/matrix.h"
#include "parsec/private_mempool.h"
#include "common.h"
#include "hicma.h"
#include "HiCMA_dpotrf_L_3flow2.h"
#include "hicma_parsec_internal.h"

static int DEBUG_INFO = 0;
double *gather_time;
double *gather_time_tmp;

/**
 *******************************************************************************
 *
 * @ingroup dplasma_complex64
 *
 * dplasma_dpotrf_New - Generates the taskpool that Computes the Cholesky
 * factorization of a symmetric positive definite (or Hermitian positive
 * definite in the complex case) matrix A, with or without recursive calls.
 * The factorization has the form
 *
 *    \f[ A = \{_{L\times L^H, if uplo = PlasmaLower}^{U^H\times U, if uplo = PlasmaUpper} \f]
 *
 *  where U is an upper triangular matrix and L is a lower triangular matrix.
 *
 * WARNING: The computations are not done by this call.
 *
 * Hierarchical DAG Scheduling for Hybrid Distributed Systems; Wu, Wei and
 * Bouteiller, Aurelien and Bosilca, George and Faverge, Mathieu and Dongarra,
 * Jack. 29th IEEE International Parallel & Distributed Processing Symposium,
 * May 2015. (https://hal.inria.fr/hal-0107835)
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          = PlasmaUpper: Upper triangle of A is referenced;
 *          = PlasmaLower: Lower triangle of A is referenced.
 *
 * @param[in,out] A
 *          Descriptor of the distributed matrix A.
 *          On exit, the uplo part of A is overwritten with the factorized
 *          matrix.
 *
 * @param[out] info
 *          Address where to store the output information of the factorization,
 *          this is not synchronized between the nodes, and might not be set
 *          when function exists.
 *          On DAG completion:
 *              - info = 0 on all nodes if successful.
 *              - info > 0 if the leading minor of order i of A is not positive
 *                definite, so the factorization could not be completed, and the
 *                solution has not been computed. Info will be equal to i on the
 *                node that owns the diagonal element (i,i), and 0 on all other
 *                nodes.
 *
 *******************************************************************************
 *
 * @return
 *          \retval NULL if incorrect parameters are given.
 *          \retval The parsec taskpool describing the operation that can be
 *          enqueued in the runtime with parsec_context_add_taskpool(). It, then, needs to be
 *          destroy with dplasma_dpotrf_Destruct();
 *
 *******************************************************************************
 *
 * @sa dplasma_dpotrf
 * @sa dplasma_dpotrf_Destruct
 * @sa dplasma_cpotrf_New
 * @sa dplasma_dpotrf_New
 * @sa dplasma_spotrf_New
 *
 ******************************************************************************/
static int wrap_potrf(parsec_execution_stream_t * es, 
                      __parsec_HiCMA_dpotrf_L_3flow2_potrf_dpotrf_task_t * this_task)
{
    parsec_HiCMA_dpotrf_L_3flow2_taskpool_t *parsec_tp = (parsec_HiCMA_dpotrf_L_3flow2_taskpool_t*)this_task->taskpool;
    /* Record start time of potrf */
    parsec_tp->_g_potrf_time_temp[0] = MPI_Wtime();
    gather_time_tmp[es->th_id] = parsec_tp->_g_potrf_time_temp[0];
    return parsec_tp->_g_wrap_potrf(es, (parsec_task_t *)this_task);
}

static int wrap_potrf_complete(parsec_execution_stream_t * es,
                      __parsec_HiCMA_dpotrf_L_3flow2_potrf_dpotrf_task_t * this_task)
{
    int val;
    double end_time;
    parsec_HiCMA_dpotrf_L_3flow2_taskpool_t *parsec_tp = (parsec_HiCMA_dpotrf_L_3flow2_taskpool_t*)this_task->taskpool;
    val = parsec_tp->_g_wrap_potrf_complete(es, (parsec_task_t *)this_task);
    end_time = MPI_Wtime();
    parsec_tp->_g_potrf_time[0] += end_time - parsec_tp->_g_potrf_time_temp[0]; 
    gather_time[es->th_id] += end_time - gather_time_tmp[es->th_id];

    if( DEBUG_INFO )
        fprintf(stderr, "POTRF %d : end_time, %lf start_time, %lf exe_time, %lf sum_time, %lf\n",
                        this_task->locals.k.value, end_time, parsec_tp->_g_potrf_time_temp[0],
                        end_time - parsec_tp->_g_potrf_time_temp[0], parsec_tp->_g_potrf_time[0]);
    return val;
}

static int wrap_trsm(parsec_execution_stream_t * es,
                      __parsec_HiCMA_dpotrf_L_3flow2_potrf_dtrsm_task_t * this_task)
{
    parsec_HiCMA_dpotrf_L_3flow2_taskpool_t *parsec_tp = (parsec_HiCMA_dpotrf_L_3flow2_taskpool_t*)this_task->taskpool;
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
                      __parsec_HiCMA_dpotrf_L_3flow2_potrf_dtrsm_task_t * this_task)
{
    parsec_HiCMA_dpotrf_L_3flow2_taskpool_t *parsec_tp = (parsec_HiCMA_dpotrf_L_3flow2_taskpool_t*)this_task->taskpool;
    int val;
    double end_time;

    if(this_task->locals.m.value == this_task->locals.k.value + 1){
        val = parsec_tp->_g_wrap_trsm_complete(es, (parsec_task_t *)this_task);
        end_time = MPI_Wtime();
        parsec_tp->_g_trsm_time[0] += end_time - parsec_tp->_g_trsm_time_temp[0];
        if( DEBUG_INFO )
            fprintf(stderr, "TRSM %d : end_time, %lf start_time, %lf exe_time, %lf sum_time, %lf\n", this_task->locals.k.value, end_time, parsec_tp->_g_trsm_time_temp[0], end_time - parsec_tp->_g_trsm_time_temp[0], parsec_tp->_g_trsm_time[0]);
    } else {
        val = parsec_tp->_g_wrap_trsm_complete(es, (parsec_task_t *)this_task);
        end_time = MPI_Wtime();
    }

    gather_time[es->th_id] += end_time - gather_time_tmp[es->th_id];
    return val;
}

static int wrap_syrk(parsec_execution_stream_t * es,
                      __parsec_HiCMA_dpotrf_L_3flow2_potrf_dsyrk_task_t * this_task)
{
    parsec_HiCMA_dpotrf_L_3flow2_taskpool_t *parsec_tp = (parsec_HiCMA_dpotrf_L_3flow2_taskpool_t*)this_task->taskpool;
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
                      __parsec_HiCMA_dpotrf_L_3flow2_potrf_dsyrk_task_t * this_task)
{
    parsec_HiCMA_dpotrf_L_3flow2_taskpool_t *parsec_tp = (parsec_HiCMA_dpotrf_L_3flow2_taskpool_t*)this_task->taskpool;
    int val;
    double end_time;

    if(this_task->locals.m.value == this_task->locals.k.value + 1){
        val = parsec_tp->_g_wrap_syrk_complete(es, (parsec_task_t *)this_task);
        end_time = MPI_Wtime();
        parsec_tp->_g_syrk_time[0] += end_time - parsec_tp->_g_syrk_time_temp[0];
        if( DEBUG_INFO )
            fprintf(stderr, "SYRK %d : end_time, %lf start_time, %lf exe_time, %lf sum_time, %lf\n", this_task->locals.k.value, end_time, parsec_tp->_g_syrk_time_temp[0], end_time - parsec_tp->_g_syrk_time_temp[0], parsec_tp->_g_syrk_time[0]);
    } else {
        val = parsec_tp->_g_wrap_syrk_complete(es, (parsec_task_t *)this_task);
        end_time = MPI_Wtime();
    }

    gather_time[es->th_id] += end_time - gather_time_tmp[es->th_id];
    return val;
}

static int wrap_gemm(parsec_execution_stream_t * es,
                      __parsec_HiCMA_dpotrf_L_3flow2_potrf_dgemm_task_t * this_task)
{
    parsec_HiCMA_dpotrf_L_3flow2_taskpool_t *parsec_tp = (parsec_HiCMA_dpotrf_L_3flow2_taskpool_t*)this_task->taskpool;
    /* Record start time of gemm */
    gather_time_tmp[es->th_id] = MPI_Wtime();
    return parsec_tp->_g_wrap_gemm(es, (parsec_task_t *)this_task);
}

static int wrap_gemm_complete(parsec_execution_stream_t * es,
                      __parsec_HiCMA_dpotrf_L_3flow2_potrf_dgemm_task_t * this_task)
{
    parsec_HiCMA_dpotrf_L_3flow2_taskpool_t *parsec_tp = (parsec_HiCMA_dpotrf_L_3flow2_taskpool_t*)this_task->taskpool;
    int val;
    double start_time = gather_time_tmp[es->th_id];
    val = parsec_tp->_g_wrap_gemm_complete(es, (parsec_task_t *)this_task);
    double end_time = MPI_Wtime();
    gather_time[es->th_id] += end_time - start_time; 
    if( DEBUG_INFO )
        fprintf(stderr, "GEMM %d %d %d : end_time, %lf start_time, %lf exe_time, %lf sum_time, %lf\n",
                        this_task->locals.m.value, this_task->locals.n.value, this_task->locals.k.value,
                        end_time, start_time, end_time - start_time, gather_time[es->th_id]);
    return val;
}

#if defined(PARSEC_HAVE_CUDA)

/* Select GPU trsm kernel 
 * Can not pass internal_taskpool, so band_size_local instead
 */
static float evaluate_gpu_potrf(parsec_task_t* task) {
    return PARSEC_HOOK_RETURN_DONE;
}

static float evaluate_gpu_trsm(parsec_task_t* task) {
    int m = ((__parsec_HiCMA_dpotrf_L_3flow2_potrf_dtrsm_task_t *)task)->locals.m.value;
    int k = ((__parsec_HiCMA_dpotrf_L_3flow2_potrf_dtrsm_task_t *)task)->locals.k.value;
    int band_size = ((__parsec_HiCMA_dpotrf_L_3flow2_potrf_dtrsm_task_t *)task)->locals.band_size_local.value;
    if( m-k < band_size )
        return PARSEC_HOOK_RETURN_DONE;
    else
        return PARSEC_HOOK_RETURN_NEXT;
}

/* Select GPU syrk kernel 
 * Can not pass internal_taskpool, so band_size_local instead
 */
static float evaluate_gpu_syrk(parsec_task_t* task) {
    int m = ((__parsec_HiCMA_dpotrf_L_3flow2_potrf_dsyrk_task_t *)task)->locals.m.value;
    int k = ((__parsec_HiCMA_dpotrf_L_3flow2_potrf_dsyrk_task_t *)task)->locals.k.value;
    int band_size = ((__parsec_HiCMA_dpotrf_L_3flow2_potrf_dsyrk_task_t *)task)->locals.band_size_local.value;
    if( 1 || m-k < band_size )
        return PARSEC_HOOK_RETURN_DONE;
    else
        return PARSEC_HOOK_RETURN_NEXT;
}

/* Select GPU gemm kernel 
 * Can not pass internal_taskpool, so band_size_local instead
 */
static float evaluate_gpu_gemm(parsec_task_t* task) {
    int m = ((__parsec_HiCMA_dpotrf_L_3flow2_potrf_dgemm_task_t *)task)->locals.m.value;
    int n = ((__parsec_HiCMA_dpotrf_L_3flow2_potrf_dgemm_task_t *)task)->locals.n.value;
    int k = ((__parsec_HiCMA_dpotrf_L_3flow2_potrf_dgemm_task_t *)task)->locals.k.value;
    int band_size = ((__parsec_HiCMA_dpotrf_L_3flow2_potrf_dgemm_task_t *)task)->locals.band_size_local.value;
    int send_full = ((__parsec_HiCMA_dpotrf_L_3flow2_potrf_dgemm_task_t *)task)->locals.send_full_tile_local.value;
    //if( m-k < band_size || (m-n < band_size && send_full) )
    //if( m-k < band_size || (n-k < band_size && m-n < band_size) )
    if( m-n < band_size )
        return PARSEC_HOOK_RETURN_DONE;
    else
        return PARSEC_HOOK_RETURN_NEXT;
}

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

parsec_taskpool_t*
HiCMA_dpotrf_L_3flow2_New( parsec_context_t *parsec,
                  int uplo,
                  parsec_tiled_matrix_dc_t *A,
                  parsec_tiled_matrix_dc_t *Ar,
                  parsec_tiled_matrix_dc_t *RG,
                  parsec_tiled_matrix_dc_t *Rank,
                  double acc, /* accuracy threshold */
                  int rk,     /* rank threshold     */
                  int storagemaxrank,  /* size of storage    */
                  int lookahead,
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
    void** eval_gpu_potrf;
    void** eval_gpu_trsm;
    void** eval_gpu_syrk;
    void** eval_gpu_gemm;

    /* Check input arguments */
    if ((uplo != PlasmaUpper) && (uplo != PlasmaLower)) {
        dplasma_error("HiCMA_dpotrf_L_3flow2_New", "illegal value of uplo");
        return NULL /*-1*/;
    }
    if ((uplo == PlasmaUpper)) {
        dplasma_error("HiCMA_dpotrf_L_3flow2_New", "HiCMA_dpotrf_L_3flow2 does not support PlasmaUpper for now");
        return NULL /*-1*/;
    }
    if(storagemaxrank > compmaxrank) {
        dplasma_error("HiCMA_dpotrf_L_3flow2_New", "maxrank for storage larger than maxrank for buffers \
                is not meaningful");
        return NULL /*-1*/;
    }
    if(storagemaxrank > (A->mb/2)) {
        dplasma_error("HiCMA_dpotrf_L_3flow2_New", "maxrank should not be larger than half of block size");
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
        parsec_HiCMA_dpotrf_L_3flow2_taskpool_t *hicma_dpotrf =
            parsec_HiCMA_dpotrf_L_3flow2_new( uplo, A, Ar, RG, Rank, info,
                                       acc, rk, storagemaxrank, compmaxrank, lookahead);

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
#if defined(PARSEC_HAVE_CUDA)
            /* GPU wrap kernel */    
            hicma_dpotrf->_g_wrap_potrf = hicma_dpotrf->super.task_classes_array[2]->incarnations[0].hook;
            hicma_dpotrf->_g_wrap_trsm = hicma_dpotrf->super.task_classes_array[3]->incarnations[0].hook;
            hicma_dpotrf->_g_wrap_syrk = hicma_dpotrf->super.task_classes_array[4]->incarnations[0].hook;
            hicma_dpotrf->_g_wrap_gemm = hicma_dpotrf->super.task_classes_array[5]->incarnations[0].hook;
            hook = (void *)&hicma_dpotrf->super.task_classes_array[2]->incarnations[0].hook;
            *hook = &wrap_potrf;
            hook = (void *)&hicma_dpotrf->super.task_classes_array[3]->incarnations[0].hook;
            *hook = &wrap_trsm;
            hook = (void *)&hicma_dpotrf->super.task_classes_array[4]->incarnations[0].hook;
            *hook = &wrap_syrk;
            hook = (void *)&hicma_dpotrf->super.task_classes_array[5]->incarnations[0].hook;
            *hook = &wrap_gemm;

            /* GPU evaluate of chores */
            eval_gpu_potrf = (void *)&hicma_dpotrf->super.task_classes_array[2]->incarnations[0].evaluate;
            eval_gpu_trsm = (void *)&hicma_dpotrf->super.task_classes_array[3]->incarnations[0].evaluate;
            eval_gpu_syrk = (void *)&hicma_dpotrf->super.task_classes_array[4]->incarnations[0].evaluate;
            eval_gpu_gemm = (void *)&hicma_dpotrf->super.task_classes_array[5]->incarnations[0].evaluate;
            *eval_gpu_potrf = &evaluate_gpu_potrf;
            *eval_gpu_trsm = &evaluate_gpu_trsm;
            *eval_gpu_syrk = &evaluate_gpu_syrk;
            *eval_gpu_gemm = &evaluate_gpu_gemm;
#endif

          /* Recursive */
        } else if( hmb < A->mb ) { 
            hicma_dpotrf->_g_wrap_potrf = hicma_dpotrf->super.task_classes_array[2]->incarnations[1].hook;
            hicma_dpotrf->_g_wrap_trsm = hicma_dpotrf->super.task_classes_array[3]->incarnations[1].hook;
            hicma_dpotrf->_g_wrap_syrk = hicma_dpotrf->super.task_classes_array[4]->incarnations[1].hook;
            hicma_dpotrf->_g_wrap_gemm = hicma_dpotrf->super.task_classes_array[5]->incarnations[1].hook;
            hook = (void *)&hicma_dpotrf->super.task_classes_array[2]->incarnations[1].hook;
            *hook = &wrap_potrf;
            hook = (void *)&hicma_dpotrf->super.task_classes_array[3]->incarnations[1].hook;
            *hook = &wrap_trsm;
            hook = (void *)&hicma_dpotrf->super.task_classes_array[4]->incarnations[1].hook;
            *hook = &wrap_syrk;
            hook = (void *)&hicma_dpotrf->super.task_classes_array[5]->incarnations[1].hook;
            *hook = &wrap_gemm;

          /* Others */
        } else {
            hicma_dpotrf->_g_wrap_potrf = hicma_dpotrf->super.task_classes_array[2]->incarnations[2].hook;
            hicma_dpotrf->_g_wrap_trsm = hicma_dpotrf->super.task_classes_array[3]->incarnations[2].hook;
            hicma_dpotrf->_g_wrap_syrk = hicma_dpotrf->super.task_classes_array[4]->incarnations[2].hook;
            hicma_dpotrf->_g_wrap_gemm = hicma_dpotrf->super.task_classes_array[5]->incarnations[2].hook;
            hook = (void *)&hicma_dpotrf->super.task_classes_array[2]->incarnations[2].hook;
            *hook = &wrap_potrf;
            hook = (void *)&hicma_dpotrf->super.task_classes_array[3]->incarnations[2].hook;
            *hook = &wrap_trsm;
            hook = (void *)&hicma_dpotrf->super.task_classes_array[4]->incarnations[2].hook;
            *hook = &wrap_syrk;
            hook = (void *)&hicma_dpotrf->super.task_classes_array[5]->incarnations[2].hook;
            *hook = &wrap_gemm;
        }

        hicma_dpotrf->_g_wrap_potrf_complete = hicma_dpotrf->super.task_classes_array[2]->complete_execution;
        hook = (void *)&hicma_dpotrf->super.task_classes_array[2]->complete_execution;
        *hook = &wrap_potrf_complete;

        hicma_dpotrf->_g_wrap_trsm_complete = hicma_dpotrf->super.task_classes_array[3]->complete_execution;
        hook = (void *)&hicma_dpotrf->super.task_classes_array[3]->complete_execution;
        *hook = &wrap_trsm_complete;

        hicma_dpotrf->_g_wrap_syrk_complete = hicma_dpotrf->super.task_classes_array[4]->complete_execution;
        hook = (void *)&hicma_dpotrf->super.task_classes_array[4]->complete_execution;
        *hook = &wrap_syrk_complete;

        hicma_dpotrf->_g_wrap_gemm_complete = hicma_dpotrf->super.task_classes_array[5]->complete_execution;
        hook = (void *)&hicma_dpotrf->super.task_classes_array[5]->complete_execution;
        *hook = &wrap_gemm_complete;

        hicma_dpotrf->_g_p_work = (parsec_memory_pool_t*)malloc(sizeof(parsec_memory_pool_t));
        //parsec_private_memory_init( hicma_dpotrf->_g_p_work, A->nb * A->nb * sizeof(double) * 8 );
        /* size used for temporary buffers obtained from hicma/compute/pzpotrf.c line 96 */
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
        //parsec_private_memory_init( hicma_dpotrf->_g_p_work, A->nb * A->nb * sizeof(double) * 8 );

        hicma_dpotrf->_g_p_work_rr = (parsec_memory_pool_t*)malloc(sizeof(parsec_memory_pool_t));
        parsec_private_memory_init( hicma_dpotrf->_g_p_work_rr, compmaxrank * compmaxrank * sizeof(double) ); 

        hicma_dpotrf->_g_p_work_mbr = (parsec_memory_pool_t*)malloc(sizeof(parsec_memory_pool_t));
        parsec_private_memory_init( hicma_dpotrf->_g_p_work_mbr, A->mb * compmaxrank * sizeof(double) ); 

        int rsvd_oversample = 10;
        int mn = rsvd_oversample + storagemaxrank;
        if(mn > A->mb)
            mn = A->mb;
        size_t rsvd_lwork = (4*mn+7) * mn;
        if(rsvd_lwork < A->mb)
            rsvd_lwork = A->mb;
        rsvd_lwork += mn*(3*A->mb+mn+1);
        size_t rsvd_liwork = 8*mn;

        hicma_dpotrf->_g_rsvd_lwork = rsvd_lwork;

        hicma_dpotrf->_g_rsvd_work = (parsec_memory_pool_t*)malloc(sizeof(parsec_memory_pool_t));
        parsec_private_memory_init( hicma_dpotrf->_g_rsvd_work, rsvd_lwork * sizeof(double) );

        hicma_dpotrf->_g_rsvd_iwork = (parsec_memory_pool_t*)malloc(sizeof(parsec_memory_pool_t));
        parsec_private_memory_init( hicma_dpotrf->_g_rsvd_iwork, rsvd_liwork * sizeof(double) );

        hicma_dpotrf->_g_d_work = (parsec_memory_pool_t*)malloc(sizeof(parsec_memory_pool_t));
        parsec_private_memory_init( hicma_dpotrf->_g_d_work, A->mb * A->mb * sizeof(double) );

        hicma_dpotrf->_g_band_size = max(((sym_two_dim_block_cyclic_band_t *)A)->band_size, 1);

        parsec_matrix_add2arena(hicma_dpotrf->arenas[PARSEC_HiCMA_dpotrf_L_3flow2_FULL_ARENA],
                                parsec_datatype_double_t, matrix_UpperLower,
                                1, A->mb, A->mb, A->mb,
                                PARSEC_ARENA_ALIGNMENT_SSE, -1 );

        parsec_matrix_add2arena(hicma_dpotrf->arenas[PARSEC_HiCMA_dpotrf_L_3flow2_DEFAULT_ARENA],
                                parsec_datatype_double_t, matrix_UpperLower,
                                1, 1, 1, 1,
                                PARSEC_ARENA_ALIGNMENT_SSE, -1 );

        parsec_matrix_add2arena(hicma_dpotrf->arenas[PARSEC_HiCMA_dpotrf_L_3flow2_UV_ARENA],
                                parsec_datatype_double_t, matrix_UpperLower,
                                1, A->mb, storagemaxrank, A->mb,
                                PARSEC_ARENA_ALIGNMENT_SSE, -1 );

        parsec_matrix_add2arena(hicma_dpotrf->arenas[PARSEC_HiCMA_dpotrf_L_3flow2_AR_ARENA],
                                parsec_datatype_int_t, matrix_UpperLower,
                                1, 1, 1, 1,
                                PARSEC_ARENA_ALIGNMENT_SSE, -1 );
    }

    return tp;
}

/**
 *******************************************************************************
 *
 * @ingroup dplasma_complex64
 *
 * dplasma_dpotrf - Computes the Cholesky factorization of a symmetric positive
 * definite (or Hermitian positive definite in the complex case) matrix A.
 * The factorization has the form
 *
 *    \f[ A = \{_{L\times L^H, if uplo = PlasmaLower}^{U^H\times U, if uplo = PlasmaUpper} \f]
 *
 *  where U is an upper triangular matrix and L is a lower triangular matrix.
 *
 *******************************************************************************
 *
 * @param[in,out] parsec
 *          The parsec context of the application that will run the operation.
 *
 * @param[in] uplo
 *          = PlasmaUpper: Upper triangle of A is referenced;
 *          = PlasmaLower: Lower triangle of A is referenced.
 *
 * @param[in] A
 *          Descriptor of the distributed matrix A.
 *          On exit, the uplo part of A is overwritten with the factorized
 *          matrix.
 *
 *******************************************************************************
 *
 * @return
 *          \retval -i if the ith parameters is incorrect.
 *          \retval 0 on success.
 *          \retval > 0 if the leading minor of order i of A is not positive
 *          definite, so the factorization could not be completed, and the
 *          solution has not been computed. Info will be equal to i on the node
 *          that owns the diagonal element (i,i), and 0 on all other nodes.
 *
 *******************************************************************************
 *
 * @sa dplasma_dpotrf_New
 * @sa dplasma_dpotrf_Destruct
 * @sa dplasma_cpotrf
 * @sa dplasma_dpotrf
 * @sa dplasma_spotrf
 *
 ******************************************************************************/
void HiCMA_dpotrf_L_3flow2_Destruct(parsec_taskpool_t* _tp)
{
    parsec_HiCMA_dpotrf_L_3flow2_taskpool_t *tp = (parsec_HiCMA_dpotrf_L_3flow2_taskpool_t*)_tp;

#if defined(PARSEC_HAVE_CUDA)
    if( tp->_g_nb_cuda_devices > 0 ) {
        workspace_memory_free( tp->_g_ws_handle );
        workspace_memory_free( tp->_g_ws_mbr );
        workspace_memory_free( tp->_g_ws_rr );

        if( NULL != tp->_g_cuda_device_index )
            free(tp->_g_cuda_device_index);
    }
#endif

    parsec_matrix_del2arena( tp->arenas[PARSEC_HiCMA_dpotrf_L_3flow2_DEFAULT_ARENA] );
    parsec_matrix_del2arena( tp->arenas[PARSEC_HiCMA_dpotrf_L_3flow2_FULL_ARENA] );
    parsec_matrix_del2arena( tp->arenas[PARSEC_HiCMA_dpotrf_L_3flow2_UV_ARENA] );
    parsec_matrix_del2arena( tp->arenas[PARSEC_HiCMA_dpotrf_L_3flow2_AR_ARENA] );
    parsec_private_memory_fini( tp->_g_p_work );
    parsec_private_memory_fini( tp->_g_p_work_mbr );
    parsec_private_memory_fini( tp->_g_p_work_rr );
    parsec_private_memory_fini( tp->_g_rsvd_work );
    parsec_private_memory_fini( tp->_g_rsvd_iwork );
    parsec_private_memory_fini( tp->_g_d_work );

    parsec_taskpool_free(_tp);
}

int
HiCMA_dpotrf_L_3flow2( parsec_context_t *parsec,
              int uplo,
              parsec_tiled_matrix_dc_t *A,
              parsec_tiled_matrix_dc_t *Ar,
              parsec_tiled_matrix_dc_t *RG,
              parsec_tiled_matrix_dc_t *Rank,
              double acc, int rk, 
              int maxrk,
              int lookahead, int hmb,
              int compmaxrank,
              int send_full_tile,
              unsigned long* tileopcounters,
              unsigned long* opcounters,
              double *critical_path_time
              )
{
    parsec_taskpool_t *hicma_zpotrf = NULL;
    int info = 0, ginfo = 0;

    /* Only for 1 vp */
    assert( parsec->nb_vp == 1 );
    int nb_threads = parsec->virtual_processes[0]->nb_cores;

    /* Allocate memory to store execution time of each process */
    gather_time = (double *)calloc(nb_threads, sizeof(double));
    gather_time_tmp = (double *)calloc(nb_threads, sizeof(double));

    hicma_zpotrf = HiCMA_dpotrf_L_3flow2_New( parsec, uplo,
                                     A, Ar, RG, Rank,
                                     acc, rk, 
                                     maxrk,
                                     lookahead,
                                     hmb,
                                     compmaxrank,
                                     send_full_tile,
                                     tileopcounters, opcounters, critical_path_time,
                                     &info );
    if( NULL != hicma_zpotrf ) {
        parsec_context_add_taskpool( parsec, (parsec_taskpool_t*)hicma_zpotrf);
        dplasma_wait_until_completion(parsec);
        HiCMA_dpotrf_L_3flow2_Destruct( hicma_zpotrf );
        if( hmb < A->mb )
            parsec_taskpool_sync_ids(); /*recursive DAGs are not synchronous on ids */
    }

    double total_time = 0.0;
    double max_time = gather_time[0];
    double min_time = gather_time[0];
    for( int i = 0; i < nb_threads; i++) {
        total_time += gather_time[i];
        if( gather_time[i] > max_time )
            max_time = gather_time[i];
        if( gather_time[i] < min_time )
            min_time = gather_time[i];
    }
   
    /* Print execution time for each process, max and min time for threads in a process */
    if( DEBUG_INFO )
        printf("gather_time %d %lf %lf %lf\n", parsec->my_rank, total_time, max_time, min_time);

    /* Free memory */
    free(gather_time);
    free(gather_time_tmp);

    /* This covers both cases when we have not compiled with MPI, or we don't need to do the reduce */
    ginfo = info;
#if defined(PARSEC_HAVE_MPI)
    /* If we don't need to reduce, don't do it, this way we don't require MPI to be initialized */
    if( A->super.nodes > 1 )
        MPI_Allreduce( &info, &ginfo, 1, MPI_INT, MPI_MAX, 
                //*(MPI_Comm*)dplasma_pcomm //TODO FIXME Allreduce of info is definitly required. Use parsec communicator
                MPI_COMM_WORLD
                );
#endif

    return ginfo;
}

