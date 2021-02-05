/**
 * @copyright (c) 2021     King Abdullah University of Science and Technology (KAUST).
 * @copyright (c) 2021     The Universiy of Tennessee and The Universiy of Tennessee Research Foundation.
 *                         All rights reserved.
 * @version 0.1.0
 * @date 2021-01-24
 *
 **/

#include "hicma_parsec.h"

/* New distribution, which is that used in
 * IPDPS2021: Leveraging parsec runtime support to tackle challenging 3d data-sparse matrix problems
 */ 

/* New rank_of for sym two dim block cyclic band */
static uint32_t hicma_parsec_sym_twoDBC_band_rank_of(parsec_data_collection_t * desc, ...)
{
    unsigned int m, n;
    va_list ap;
    sym_two_dim_block_cyclic_band_t * dc = (sym_two_dim_block_cyclic_band_t *)desc;

    /* Get coordinates */
    va_start(ap, desc);
    m = va_arg(ap, unsigned int);
    n = va_arg(ap, unsigned int);
    va_end(ap);

    /* Check tile location within band_size */
    if( (unsigned int)abs((int)m-(int)n) < dc->band_size ) {
        /* New index */
        unsigned int m_band = (unsigned int)abs((int)m - (int)n);
        return dc->band.super.super.rank_of(&dc->band.super.super, m_band, m);
    }

    return dc->off_band.super.super.rank_of(&dc->off_band.super.super, m, n);
}

/* New vpid_of for sym two dim block cyclic band */
static int32_t hicma_parsec_sym_twoDBC_band_vpid_of(parsec_data_collection_t * desc, ...)
{
    unsigned int m, n;
    va_list ap;
    sym_two_dim_block_cyclic_band_t * dc = (sym_two_dim_block_cyclic_band_t *)desc;

    /* Get coordinates */
    va_start(ap, desc);
    m = va_arg(ap, unsigned int);
    n = va_arg(ap, unsigned int);
    va_end(ap);

    /* Check tile location within band_size */
    if( (unsigned int)abs((int)m - (int)n) < dc->band_size ) {
        /* The new m in band */
        unsigned int m_band = (unsigned int)abs((int)m - (int)n);
        return dc->band.super.super.vpid_of(&dc->band.super.super, m_band, m);
    }

    return dc->off_band.super.super.vpid_of(&dc->off_band.super.super, m, n);
}

/* New data_of for sym two dim block cyclic band */
static parsec_data_t* hicma_parsec_sym_twoDBC_band_data_of(parsec_data_collection_t *desc, ...)
{
    unsigned int m, n;
    va_list ap;
    sym_two_dim_block_cyclic_band_t * dc;
    dc = (sym_two_dim_block_cyclic_band_t *)desc;

    /* Get coordinates */
    va_start(ap, desc);
    m = va_arg(ap, unsigned int);
    n = va_arg(ap, unsigned int);
    va_end(ap);

#if defined(DISTRIBUTED)
    assert(desc->myrank == hicma_parsec_sym_twoDBC_band_rank_of(desc, m, n));
#endif

    /* Check tile location within band_size */
    if( (unsigned int)abs((int)m - (int)n) < dc->band_size ) {
        /* The new m in band */
        unsigned int m_band = (unsigned int)abs((int)m - (int)n);
        return dc->band.super.super.data_of(&dc->band.super.super, m_band, m);
    }

    return dc->off_band.super.super.data_of(&dc->off_band.super.super, m, n);
}

/* New rank_of_key for sym two dim block cyclic band */
static uint32_t hicma_parsec_sym_twoDBC_band_rank_of_key(parsec_data_collection_t *desc, parsec_data_key_t key)
{
    int m, n;
    twoDBC_key_to_coordinates(desc, key, &m, &n);
    return hicma_parsec_sym_twoDBC_band_rank_of(desc, m, n);
}

/* New vpid_of_key for two dim block cyclic band */
static int32_t hicma_parsec_sym_twoDBC_band_vpid_of_key(parsec_data_collection_t *desc, parsec_data_key_t key)
{
    int m, n;
    twoDBC_key_to_coordinates(desc, key, &m, &n);
    return hicma_parsec_sym_twoDBC_band_vpid_of(desc, m, n);
}

/* New data_of_key for sym two dim block cyclic band */
static parsec_data_t* hicma_parsec_sym_twoDBC_band_data_of_key(parsec_data_collection_t *desc, parsec_data_key_t key)
{
    int m, n;
    twoDBC_key_to_coordinates(desc, key, &m, &n);
    return hicma_parsec_sym_twoDBC_band_data_of(desc, m, n);
}

/* 
 * sysm_two_dim_block_cyclic_band_t structure init 
 * It inherits from off-band, so should be called after initialization of off_band
 */
void hicma_parsec_sym_two_dim_block_cyclic_band_init( sym_two_dim_block_cyclic_band_t *desc,
                                         int nodes, int myrank, int band_size ) {
    parsec_tiled_matrix_dc_t *off_band = &desc->off_band.super;
    parsec_data_collection_t *dc = (parsec_data_collection_t*)desc;

    parsec_tiled_matrix_dc_init( &desc->super, off_band->mtype, off_band->storage, off_band->dtype,
                                 nodes, myrank, off_band->mb, off_band->nb, off_band->lm, off_band->ln,
                                 off_band->i, off_band->j, off_band->m, off_band->n );

    desc->band_size  = band_size;
    dc->rank_of      = hicma_parsec_sym_twoDBC_band_rank_of;
    dc->vpid_of      = hicma_parsec_sym_twoDBC_band_vpid_of;
    dc->data_of      = hicma_parsec_sym_twoDBC_band_data_of;
    dc->rank_of_key  = hicma_parsec_sym_twoDBC_band_rank_of_key;
    dc->vpid_of_key  = hicma_parsec_sym_twoDBC_band_vpid_of_key;
    dc->data_of_key  = hicma_parsec_sym_twoDBC_band_data_of_key;
}
