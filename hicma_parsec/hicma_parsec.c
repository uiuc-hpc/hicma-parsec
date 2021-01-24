/**
 * @copyright (c) 2021 King Abdullah University of Science and Technology (KAUST).
 *                     The Universiy of Tennessee and The Universiy of Tennessee Research Foundation.
 *                     All rights reserved.
 **/

#include "hicma_parsec.h"

#ifdef PARSEC_HAVE_UNISTD_H
#include <unistd.h>
#endif
#ifdef PARSEC_HAVE_LIMITS_H
#include <limits.h>
#endif
#if defined(PARSEC_HAVE_GETOPT_H)
#include <getopt.h>
#endif  /* defined(PARSEC_HAVE_GETOPT_H) */

/*******************************
 * globals and argv set values *
 *******************************/
/* Gather time */
double time_elapsed = 0.0;
double sync_time_elapsed = 0.0;
double *gather_time;
double *gather_time_tmp;

#if DYNAMIC_COLLECTIVE_PATTERN
/* Parsec routines for collective pattern */
extern int remote_dep_bcast_star_child(int me, int him);
extern int (*remote_dep_bcast_child)(int me, int him);
#endif

/**********************************
 * Command line arguments
 **********************************/
static void print_usage(void)
{
	fprintf(stderr,
			"Mandatory argument:\n"
			" number            : dimension (N) of the matrices (required)\n"
			"Optional arguments:\n"
			" -P --grid-rows : rows (P) in the PxQ process grid   (default: NP)\n"
			" -Q --grid-cols : columns (Q) in the PxQ process grid (default: NP/P)\n"
			"\n"
			" -N                : dimension (N) of the matrices (required)\n"
			" -t --NB           : tile size  (required)\n"
			" -z --HNB          : Inner NB used for recursive algorithms (default: MB)\n"
			" -x --check        : verify the results\n"
			" -v --verbose      : extra verbose output\n"
			" -c --cores        : number of concurent threads (default: number of physical hyper-threads)\n"
			" -g --gpus         : number of GPU (default: 0)\n"
			" -f --fixedrank        : Fixed rank threshold used in recompression stage of HCORE_GEMM\n"
			" -e --fixedacc         : Fixed accuracy threshold used in recompression stage of HCORE_GEMM\n"
			"\n"
			" -u --maxrank          : Maxrank limit used in creation of descriptors\n"
			" -G --genmaxrank       : Maxrank limit used in generation\n"
			" -U --compmaxrank      : Maxrank limit used in allocation of buffers for HiCMA_dpotrf operation.\n"
			"NOTE: These three maxrank values might be equal."
			"\n"
			" -w --wavek            : Wave number for electrodynamics problem\n"
			" -j --adddiag          : Add this number to diagonal elements to make the matrix positive definite in electrodynamics problem\n"
			" -Z --band             : if 0, normal two_dim_block_cyclic, but not contiguous memory allocation for whole matrix; if > 0, with band distribution \n"
			" -Y --lookahead        : set lookahead, from -1 to NT-1; default -1, will set to auto_tuned band_size \n"
			" -E --auto-band        : set auto select the most suitable band size \n"
			" -D --kind_of_problem  : 0: randtlr\n"
			"                         1: electrodynamics-2d-sinus\n"
			"                         2: statistics-2d-sqexp \n"
			"                         3: statistics-3d-sqexp \n"
            "                         4: statistics-3d-exp \n"
            "                         5: electrodynamics-3d-sin \n"
            "                         6: md-3d-virus \n"
            "                         7: md-3d-cube \n"
            " -M --mesh_file       : path to mesh file\n"
            " -K --numobj          : numnber of objects (number of viurses within a population)\n"
            " -B --rbf_kernel      : type of RBF basis function (0:Gaussian, 1:Expon, 2:InvQUAD, 3:InvMQUAD, 4:Maternc1, 5:Maternc2, 6:TPS, 7:CTPS, and 8:Wendland)\n"
            " -R --radius          : radius of influential nodes\n"
            " -O --order           : No, Morton, or Hilbert ordering (0, 1, or 2, respectively) starsh supports only Morton ordering. For Hilbert ordering, you need to use python code by importing Hilbert package. SOON will be provides by starsh\n"
            " -S --density               : denstiy of sphere packing for rbf application. if you want random distibution, type -1 \n"
            " -h --help            : this message\n"
			"\n");
	parsec_usage();
}

#define GETOPT_STRING "c:P:Q:N:t:h::f:e:u:G:U::w:j:Z:D:F:Y:E:W:M:K:B:R:O:S"

#if defined(PARSEC_HAVE_GETOPT_LONG)
static struct option long_options[] =
{
	/* PaRSEC specific options */
	{"cores",       required_argument,  0, 'c'},
	{"c",           required_argument,  0, 'c'},

	/* Generic Options */
	{"grid-rows",   required_argument,  0, 'P'},
	{"P",           required_argument,  0, 'P'},
	{"grid-cols",   required_argument,  0, 'Q'},
	{"Q",           required_argument,  0, 'Q'},

	{"N",           required_argument,  0, 'N'},
	{"NB",          required_argument,  0, 't'},
	{"t",           required_argument,  0, 't'},
	{"check",       no_argument,        0, 'x'},
	{"x",           no_argument,        0, 'x'},

	/* Recursive options */
	{"z",           required_argument,  0, 'z'},
	{"HNB",         required_argument,  0, 'z'},

	/* Auxiliary options */
	{"verbose",     optional_argument,  0, 'v'},
	{"v",           optional_argument,  0, 'v'},
	{"help",        no_argument,        0, 'h'},
	{"h",           no_argument,        0, 'h'},

	/* HiCMA options */
	{"fixedrank",   required_argument,  0, 'f'},
	{"f",           required_argument,  0, 'f'},
	{"fixedacc",    required_argument,  0, 'e'},
	{"e",           required_argument,  0, 'e'},
	{"maxrank",     required_argument,  0, 'u'},
	{"u",           required_argument,  0, 'u'},
	{"genmaxrank",  required_argument,  0, 'G'},
	{"G",           required_argument,  0, 'G'},
	{"compmaxrank", required_argument,  0, 'U'},
	{"U",           required_argument,  0, 'U'},
	{"wavek",       required_argument,  0, 'w'},
	{"w",           required_argument,  0, 'w'},
	{"adddiag",     required_argument,  0, 'j'},
	{"j",           required_argument,  0, 'j'},
	{"band",        required_argument,  0, 'Z'},
	{"Z",           required_argument,  0, 'Z'},
	{"lookahead",   required_argument,  0, 'Y'},
	{"Y",           required_argument,  0, 'Y'},
	{"kind_of_problem",        required_argument,  0, 'D'},
	{"D",           required_argument,  0, 'D'},
	{"send_full_tile",        required_argument,  0, 'F'},
	{"F",           required_argument,  0, 'F'},
	{"auto-band",   required_argument,  0, 'E'},
	{"E",           required_argument,  0, 'E'},
	{"two-flow",    required_argument,  0, 'W'},
    {"W",           required_argument,  0, 'W'},
    {"mesh_file",    required_argument, 0, 'M'},
    {"M",            required_argument,  0, 'M'},
    {"numobj",       required_argument, 0, 'K'},
    {"K",           required_argument,  0, 'K'},
    {"rbf_kernel",  required_argument, 0, 'B'},
    {"B",           required_argument,  0, 'B'},
    {"radius",     required_argument, 0, 'R'},
    {"R",           required_argument,  0, 'R'},
    {"order",    required_argument, 0, 'O'},
    {"O",           required_argument,  0, 'O'},
    {"density",    required_argument, 0, 'S'},
    {"S",           required_argument,  0, 'S'},
    {0, 0, 0, 0}
};
#endif  /* defined(PARSEC_HAVE_GETOPT_LONG) */

static void parse_arguments(int *_argc, char*** _argv, int* iparam, double* dparam)
{
    extern char **environ;
    int opt = 0;
    int rc, c;
	int argc = *_argc;
	char **argv = *_argv;
	char *add_dot = NULL;
	char *value;

    /* Default */
    iparam[IPARAM_NCORES] = -1;         // number_of_cores_per_node - 1
    iparam[IPARAM_NGPUS]  = 0;          // no GPU
    iparam[IPARAM_BAND] = 1;            // band_size 1
    iparam[IPARAM_LOOKAHEAD] = -1;      // lookahead set to auto-tuned band_size
    iparam[IPARAM_KIND_OF_PROBLEM] = 2; // statistics-2d-sqexp 
    iparam[IPARAM_SEND_FULL_TILE]  = 0; // do not send full tile
    iparam[IPARAM_HNB] = 300;           // subtile size in recursive
    iparam[IPARAM_AUTO_BAND] = 1;       // enable auto-tuning band_size
    iparam[IPARAM_TWO_FLOW] = 0;        // not force to run 2flow version
    iparam[IPARAM_P] = 0;               // row process grid, with set to numberi_of_nodes if not set
    iparam[IPARAM_N] = 0;               // matrix size, need to set and will check later 
    iparam[IPARAM_NB] = 0;              // tile size, need to set and will check later 
    iparam[IPARAM_VERBOSE] = 0;         // verbose
    dparam[DPARAM_ADD_DIAG] = 0.0;      // set to matrix size
    dparam[DPARAM_WAVEK] = 50;          // wave_k in synthetic 2D application (-D 1)
    dparam[DPARAM_FIXED_ACC] = 1.0e-8;  // default yield 1.0e-9
    iparam[IPARAM_MAX_RANK] = 0;        // maxrank, set to tile_size/2 by default
    iparam[IPARAM_GEN_MAX_RANK] = 0;    // default IPARAM_MAX_RANK
    iparam[IPARAM_COMP_MAX_RANK] = 0;   // default IPARAM_MAX_RANK


	do {
#if defined(PARSEC_HAVE_GETOPT_LONG)
		c = getopt_long_only(argc, argv, "",
				long_options, &opt);
#else
		c = getopt(argc, argv, GETOPT_STRING);
		(void) opt;
#endif  /* defined(PARSEC_HAVE_GETOPT_LONG) */

		// printf("%c: %s = %s\n", c, long_options[opt].name, optarg);
		switch(c)
		{
			case 'c': iparam[IPARAM_NCORES] = atoi(optarg); break;
			case 'P': iparam[IPARAM_P] = atoi(optarg); break;
			case 'Q': iparam[IPARAM_Q] = atoi(optarg); break;
			case 'N': iparam[IPARAM_N] = atoi(optarg); break;
			case 't': iparam[IPARAM_NB] = atoi(optarg); break;
			case 'x': iparam[IPARAM_CHECK] = 1; iparam[IPARAM_VERBOSE] = 1; break;
			case 'z': iparam[IPARAM_HNB] = atoi(optarg); break;
            case 'f': iparam[IPARAM_FIXED_RANK]  = atoi(optarg); break;
            case 'u': iparam[IPARAM_MAX_RANK]  = atoi(optarg); break;
            case 'G': iparam[IPARAM_GEN_MAX_RANK]  = atoi(optarg); break;
            case 'U': iparam[IPARAM_COMP_MAX_RANK]  = atoi(optarg); break;
            case 'w': dparam[DPARAM_WAVEK]  = atof(optarg); break;
            case 'e': dparam[DPARAM_FIXED_ACC]  = atof(optarg); break;
			case 'j': dparam[DPARAM_ADD_DIAG]  = atof(optarg); break;
			case 'Z': iparam[IPARAM_BAND]  = atoi(optarg); break;
			case 'Y': iparam[IPARAM_LOOKAHEAD]  = atoi(optarg); break;
			case 'D': iparam[IPARAM_KIND_OF_PROBLEM]  = atoi(optarg); break;
			case 'F': iparam[IPARAM_SEND_FULL_TILE]  = atoi(optarg); break;
			case 'E': iparam[IPARAM_AUTO_BAND]  = atoi(optarg); break;
            case 'W': iparam[IPARAM_TWO_FLOW]  = atoi(optarg); break;
            case 'M' : mesh_file = optarg; break;
            case 'K': iparam[IPARAM_NUMOBJ]  = atoi(optarg); break;
            case 'B': iparam[IPARAM_RBFKERNEL]  = atoi(optarg); break;
            case 'O': iparam[IPARAM_ORDER]  = atoi(optarg); break;
            case 'R': dparam[DPARAM_RAD]    = atof(optarg); break;
            case 'S': dparam[DPARAM_DENST]  = atof(optarg); break;
            case 'h': print_usage(); exit(0); break;

            case 'v':
                      if(optarg)  iparam[IPARAM_VERBOSE] = atoi(optarg);
                      else        iparam[IPARAM_VERBOSE] = 1;
                      break;

            case '?': /* getopt_long already printed an error message. */
                      exit(1);
                      break;  /* because some compilers are just too annoying */

            default:
                      break; /* Assume anything else is parsec/mpi stuff */
        }
    } while(-1 != c);

    int verbose = iparam[IPARAM_RANK] ? 0 : iparam[IPARAM_VERBOSE];

	if(iparam[IPARAM_NGPUS] < 0) iparam[IPARAM_NGPUS] = 0;
	if(iparam[IPARAM_NGPUS] > 0) {
		if( iparam[IPARAM_VERBOSE] ) {
			parsec_setenv_mca_param( "device_show_capabilities", "1", &environ );
			parsec_setenv_mca_param( "device_show_statistics", "1", &environ );
		}
	}

	/* Check the process grid */
	if(0 == iparam[IPARAM_P])
		iparam[IPARAM_P] = iparam[IPARAM_NNODES];
	else if(iparam[IPARAM_P] > iparam[IPARAM_NNODES])
	{
		fprintf(stderr, "#XXXXX There are only %d nodes in the world, and you requested P=%d\n",
				iparam[IPARAM_NNODES], iparam[IPARAM_P]);
		exit(2);
	}
	if(0 == iparam[IPARAM_Q])
		iparam[IPARAM_Q] = iparam[IPARAM_NNODES] / iparam[IPARAM_P];
	int pqnp = iparam[IPARAM_Q] * iparam[IPARAM_P];
	if(pqnp > iparam[IPARAM_NNODES])
	{
		fprintf(stderr, "#XXXXX the process grid PxQ (%dx%d) is larger than the number of nodes (%d)!\n", iparam[IPARAM_P], iparam[IPARAM_Q], iparam[IPARAM_NNODES]);
		exit(2);
	}
	if(verbose && (pqnp < iparam[IPARAM_NNODES]))
	{
		fprintf(stderr, "#!!!!! the process grid PxQ (%dx%d) is smaller than the number of nodes (%d). Some nodes are idling!\n", iparam[IPARAM_P], iparam[IPARAM_Q], iparam[IPARAM_NNODES]);
	}

	/* Matrix size is reqired !! */ 
	if( iparam[IPARAM_N] <= 0 )
		fprintf(stderr, "#XXXXX the matrix size (N) is not set!\n");
                
	/* Tile size is reqired !! */ 
	if( iparam[IPARAM_NB] <= 0 )
		fprintf(stderr, "#XXXXX the tile size (NB) is not set!\n");

        /* Set add_diag to matrix_size */
        if( dparam[DPARAM_ADD_DIAG] <= 0.0 )
            dparam[DPARAM_ADD_DIAG] = (double)iparam[IPARAM_N];

        /* If maxrank not set */
	if( iparam[IPARAM_MAX_RANK] <= 0 ) {
		iparam[IPARAM_MAX_RANK] = iparam[IPARAM_NB] / 2;
		fprintf(stderr, "Max rank has not been specified. Forced to the size of a tile %d\n",
				iparam[IPARAM_MAX_RANK]);
	}

        /* Set generated maxrank to maxrank by default */
        if( iparam[IPARAM_GEN_MAX_RANK] <= 0 )
		iparam[IPARAM_GEN_MAX_RANK] = iparam[IPARAM_MAX_RANK];

        /* Set computed maxrank to maxrank by default */
        if( iparam[IPARAM_COMP_MAX_RANK] <= 0 )
		iparam[IPARAM_COMP_MAX_RANK] = iparam[IPARAM_MAX_RANK];

	(void)rc;
}

static void print_arguments(int* iparam)
{
	int verbose = iparam[IPARAM_RANK] ? 0 : iparam[IPARAM_VERBOSE];

	if( verbose ) {
		fprintf(stderr, "#+++++ cores detected       : %d\n", iparam[IPARAM_NCORES]);

		fprintf(stderr, "#+++++ nodes x cores + gpu  : %d x %d + %d (%d+%d)\n"
				"#+++++ P x Q                : %d x %d (%d/%d)\n",
				iparam[IPARAM_NNODES],
				iparam[IPARAM_NCORES],
				iparam[IPARAM_NGPUS],
				iparam[IPARAM_NNODES] * iparam[IPARAM_NCORES],
				iparam[IPARAM_NNODES] * iparam[IPARAM_NGPUS],
				iparam[IPARAM_P], iparam[IPARAM_Q],
				iparam[IPARAM_Q] * iparam[IPARAM_P], iparam[IPARAM_NNODES]);

		fprintf(stderr, "#+++++ M x N                : %d x %d\n",
				iparam[IPARAM_N], iparam[IPARAM_N]);

		fprintf(stderr, "#+++++ MB x NB              : %d x %d\n",
				iparam[IPARAM_NB], iparam[IPARAM_NB]);
		fprintf(stderr, "#+++++ HMB x HNB            : %d x %d\n", iparam[IPARAM_HNB], iparam[IPARAM_HNB]);
	}
}

parsec_context_t* setup_parsec(int argc, char **argv, int *iparam, double *dparam)
{
#ifdef PARSEC_HAVE_MPI
	{
		int provided;
		MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
	}
	MPI_Comm_size(MPI_COMM_WORLD, &iparam[IPARAM_NNODES]);
	MPI_Comm_rank(MPI_COMM_WORLD, &iparam[IPARAM_RANK]);
#else
	iparam[IPARAM_NNODES] = 1;
	iparam[IPARAM_RANK] = 0;
#endif
	parse_arguments(&argc, &argv, iparam, dparam);
	int verbose = iparam[IPARAM_VERBOSE];
	if(iparam[IPARAM_RANK] > 0) verbose = 0;

	/* Timer start */
	TIME_START();

	/* Once we got out arguments, we should pass whatever is left down */
	int parsec_argc, idx;
	char** parsec_argv = (char**)calloc(argc, sizeof(char*));
	parsec_argv[0] = argv[0];  /* the app name */
	for( idx = parsec_argc = 1;
			(idx < argc) && (0 != strcmp(argv[idx], "--")); idx++);
	if( idx != argc ) {
		for( parsec_argc = 1, idx++; idx < argc;
				parsec_argv[parsec_argc] = argv[idx], parsec_argc++, idx++);
	}
	parsec_context_t* ctx = parsec_init(iparam[IPARAM_NCORES],
			&parsec_argc, &parsec_argv);
	free(parsec_argv);
	if( NULL == ctx ) {
		/* Failed to correctly initialize. In a correct scenario report
		 * upstream, but in this particular case bail out.
		 */
		exit(-1);
	}

	/* If the number of cores has not been defined as a parameter earlier
	   update it with the default parameter computed in parsec_init. */
	if(iparam[IPARAM_NCORES] <= 0)
	{
		int p, nb_total_comp_threads = 0;
		for(p = 0; p < ctx->nb_vp; p++) {
			nb_total_comp_threads += ctx->virtual_processes[p]->nb_cores;
		}
		iparam[IPARAM_NCORES] = nb_total_comp_threads;
	}
	print_arguments(iparam);

        /* Check N and NB setted */
        if( 0 == iparam[IPARAM_N] || 0 == iparam[IPARAM_NB] )
            exit(1);

	if( verbose ) TIME_PRINT(iparam[IPARAM_RANK], ("PaRSEC initialized\n"));
	return ctx;
}

void cleanup_parsec(parsec_context_t* parsec, int *iparam)
{
	parsec_fini(&parsec);

#ifdef PARSEC_HAVE_MPI
	MPI_Finalize();
#endif
	(void)iparam;
}

/* Main routine for TLR Cholesky
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
 * @param [in] two_flow:           force to run two_flow version even if band_size == 1 
 * @param [in] tileopcounters:     count the number of tiles 
 * @param [in] opcounters:         count operations during factorization 
 * @param [in] critical_path_time: measure the critical path time 
 * @return info:                   0 on all nodes if successful. 
 *                                 > 0 if the leading minor of order i of A is not positive
 *                                 definite, so the factorization could not be completed, and the
 *                                 solution has not been computed. Info will be equal to i on the
 *                                 node that owns the diagonal element (i,i), and 0 on all other nodes
 */
int HiCMA_dpotrf_L( parsec_context_t *parsec,
		int uplo,
		parsec_tiled_matrix_dc_t *A,
		parsec_tiled_matrix_dc_t *Ar,
		parsec_tiled_matrix_dc_t *Rank,
		double acc, int rk,
		int maxrank,
		int *lookahead,
		int band_size,
		int hmb,
		int compmaxrank,
		int send_full_tile,
		int *two_flow, 
		unsigned long* tileopcounters,
		unsigned long* opcounters,
		double *critical_path_time,
		int verbose
		)
{
	parsec_taskpool_t *hicma_zpotrf = NULL;
	int info = 0, ginfo = 0;

	/* Only for 1 vp */
	assert( parsec->nb_vp == 1 );
	int nb_threads = parsec->virtual_processes[0]->nb_cores;

        /* Set lookahead to auto-tuned band_size if lookahead = -1 */
        assert( *lookahead >= -1 );
        if( -1 == *lookahead ) {
            *lookahead = band_size;
	    if( 0 == A->super.myrank && verbose )
                printf("Lookahead is not provided, set lookahead = band_size = %d\n", *lookahead);
        }

#if !DYNAMIC_COLLECTIVE_PATTERN
        /* TIPS for performance */
 	if( 0 == A->super.myrank && band_size > 1 ) { 
 		fprintf(stderr, YEL "\nWARNING: band_size= %d (> 1), so add flag '-- -mca runtime_comm_coll_bcast 0' at the end of command for better performance !!!\n\n" RESET, band_size);
 	}
#endif

	/* Allocate memory to store execution time of each process */
	gather_time = (double *)calloc(nb_threads, sizeof(double));
	gather_time_tmp = (double *)calloc(nb_threads, sizeof(double));

	/* Call 3flow version
	 *
	 * This version is used in papers of ProTools at SC2019 and PASC2020, 
	 * only diagonal is dense, plus dynamical memory allocation based on the actual rank
	 *
	 * ProTools at SC2019: Performance analysis of tile low-rank cholesky factorization using parsec instrumentation tools
	 * PASC2020: Extreme-scale task-based cholesky factorization toward climate and weather prediction applications
	 *
	 */ 
	if( 1 == band_size && !(*two_flow) ) {
		int nodes = A->super.nodes;
		int rank = A->super.myrank;
		int MB = A->mb;
		int LDA = A->m;
		int N_UV = maxrank * A->lnt;
		int P = ((two_dim_block_cyclic_t *)Ar)->grid.rows;

		if( 0 == A->super.myrank && verbose )
			printf(MAG "3flow version start\n" RESET);

		/* dcAv contains V. */
		sym_two_dim_block_cyclic_t dcAv;
		sym_two_dim_block_cyclic_init(&dcAv, matrix_RealDouble,
				rank, MB, maxrank, LDA, N_UV, 0, 0,
				LDA, N_UV, P, nodes/P, uplo);
		parsec_data_collection_set_key((parsec_data_collection_t*)&dcAv, "dcAv");

		/* Set address in Av to A+MB*maxrank*sizeof(double) */
		parsec_tiled_matrix_dc_t *Av = (parsec_tiled_matrix_dc_t *)&dcAv;
		parsec_Av_memory( parsec, A, Av, Ar, maxrank );

		hicma_zpotrf = HiCMA_dpotrf_L_3flow_New( parsec, uplo,
				A, A, Av, Ar, Rank,
				acc, rk,
				maxrank,
				*lookahead,
				band_size,
				hmb,
				compmaxrank,
				send_full_tile,
				tileopcounters, opcounters, critical_path_time,
				&info );

		if( NULL != hicma_zpotrf ) {
			parsec_context_add_taskpool( parsec, (parsec_taskpool_t*)hicma_zpotrf);
			parsec_context_start(parsec);
			parsec_context_wait(parsec);
			HiCMA_dpotrf_L_3flow_Destruct( hicma_zpotrf );
			if( hmb < A->mb )
				parsec_taskpool_sync_ids(); /*recursive DAGs are not synchronous on ids */
		}

	} else {
		/* The 2flow version
		 *
		 * Combine U and V flow in the two flow version
		 * Hybrid dense and TLR tiles
		 *
		 * Could also be used when band_size = 1, and difference between 3flow version is the way deal with U and V 
		 *
		 * IPDPS2021: Leveraging parsec runtime support to tackle challenging 3d data-sparse matrix problems
		 *
		 */
                if( 0 == A->super.myrank && verbose )
                        printf(BLU "2flow version start\n" RESET);
		*two_flow = 1;

		/* Set the right collective pattern
		 * need to be one taskpoll at a time */
#if DYNAMIC_COLLECTIVE_PATTERN
		if( band_size > 1 ) {
			remote_dep_bcast_child = remote_dep_bcast_star_child;
		}
#endif

		hicma_zpotrf = HiCMA_dpotrf_L_2flow_New( parsec, uplo,
				A, Ar, Rank,
				acc, rk,
				maxrank,
				*lookahead,
				band_size,
				hmb,
				compmaxrank,
				send_full_tile,
				tileopcounters, opcounters, critical_path_time,
				&info );

		if( NULL != hicma_zpotrf ) {
			parsec_context_add_taskpool( parsec, (parsec_taskpool_t*)hicma_zpotrf);
			parsec_context_start(parsec);
			parsec_context_wait(parsec);
			HiCMA_dpotrf_L_2flow_Destruct( hicma_zpotrf );
			if( hmb < A->mb )
				parsec_taskpool_sync_ids(); /*recursive DAGs are not synchronous on ids */
		}
	}

#if PRINT_THREAD_EXE_TIME 
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
	fprintf(stderr, "Execution_time_each_process %d : %lf %lf %lf\n", parsec->my_rank, total_time, max_time, min_time);

	/* Print execution time for each thread */
	for( int i = 0; i < nb_threads; i++ )
		fprintf(stderr, "Execution_time_each_thread %d %d : %lf\n",
				parsec->my_rank, i, gather_time[i]);
#endif

	/* Free memory */
	free(gather_time);
	free(gather_time_tmp);

	/* This covers both cases when we have not compiled with MPI, or we don't need to do the reduce */
	ginfo = info;
#if defined(PARSEC_HAVE_MPI)
	/* If we don't need to reduce, don't do it, this way we don't require MPI to be initialized */
	if( A->super.nodes > 1 )
		MPI_Allreduce( &info, &ginfo, 1, MPI_INT, MPI_MAX,
				MPI_COMM_WORLD
			     );
#endif

	return ginfo;
}
