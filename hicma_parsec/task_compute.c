/**
 * @copyright (c) 2020 King Abdullah University of Science and Technology (KAUST).
 *                     The Universiy of Tennessee and The Universiy of Tennessee Research Foundation.
 *                     All rights reserved.
 **/
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>

/* Calculate tasks for each node 
 * Assume load balanced
 */
int main(int argc, char *argv[]){
    long long int NT, matrix_size = 16;
    int tile_size = 4, nodes = 1, ch, rank = 2, verbose = 0; 
    long double po, trsm, syrk, gemm, count;

    while ((ch = getopt(argc, argv, "N:t:u:a:v:")) != -1) {
        switch (ch) {
            case 'N': matrix_size = atoi(optarg); break;
            case 't': tile_size = atoi(optarg); break;
            case 'u': rank = atoi(optarg); break;
            case 'a': nodes = atoi(optarg); break;
            case 'v': verbose = 1; break;
            case '?': case 'h': default:
                fprintf(stderr,
                        "SUPER:\n"
                        "-N : dimension (N) of the matrices (default: 16)\n"
                        "-t : dimension (NB) of the tiles (default: 4)\n"
                        "-u : rank of the tiles (default: 2)\n"
                        "-a : No. of nodes\n"
                        "-v : more outputs\n"
                        "\n");
            exit(1);
        }
    }

    NT = (matrix_size % tile_size)? (matrix_size/tile_size + 1): matrix_size/tile_size;
    po = NT;
    trsm = NT * (NT-1) / 2;
    syrk = NT * (NT-1) / 2;
    gemm = NT * (NT-1) * (NT-2) / 6;
    count = po * tile_size * tile_size * tile_size / 3 + trsm * tile_size * tile_size * rank + syrk * 2 * (tile_size * tile_size * rank + tile_size * rank * rank * 2);
    if( verbose ) {
        printf("Matrix size %lld of tile size %d, No. of %d nodes, each node there are %llfd po\n", matrix_size, tile_size, nodes, po/nodes); 
        printf("Matrix size %lld of tile size %d, No. of %d nodes, each node there are %llf trsm\n", matrix_size, tile_size, nodes, trsm/nodes); 
        printf("Matrix size %lld of tile size %d, No. of %d nodes, each node there are %llf syrk\n", matrix_size, tile_size, nodes, syrk/nodes); 
        printf("Matrix size %lld of tile size %d, No. of %d nodes, each node there are %llf gemm\n", matrix_size, tile_size, nodes, gemm/nodes); 
    }

    printf("Matrix size %lld of tile size %d, No. of %d nodes, rank of %d: operations, total: %llg, each node %llg\n", matrix_size, tile_size, nodes, rank, count, count/nodes); 

    return 0;
}
