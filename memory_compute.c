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

#define MEM 128

/* Calculate memory/nodes needed */
int main(int argc, char *argv[]){
    long long int memory, memory_nodes, NT, matrix_size = 16;
    int tile_size = 4, maxrank = 4, nodes = 0, ch; 
    double uv, Ag;

    while ((ch = getopt(argc, argv, "N:t:u:a:")) != -1) {
        switch (ch) {
            case 'N': matrix_size = atoi(optarg); break;
            case 't': tile_size = atoi(optarg); break;
            case 'u': maxrank = atoi(optarg); break;
            case 'a': nodes = atoi(optarg); break;
            case '?': case 'h': default:
                fprintf(stderr,
                        "SUPER:\n"
                        "-N : dimension (N) of the matrices (default: 16)\n"
                        "-t : dimension (NB) of the tiles (default: 4)\n"
                        "-u : maxrank of the tiles (default: 4)\n"
                        "-a : No. of nodes\n"
                        "\n");
            exit(1);
        }
    }

    memory = MEM * 1024ll * 1024 * 1024;
    memory_nodes = nodes * 1024ll * 1024 * 1024;
    NT = (matrix_size % tile_size)? (matrix_size/tile_size + 1): matrix_size/tile_size;
    Ag = (double)NT * tile_size * tile_size * 8;
    uv = ((double)NT * (NT-1)) * tile_size * maxrank * 8; 

    if( nodes )
        printf("Matrix size %lld of tile size %d, No. of %d nodes, if maxrank is %d each node requries %lf G\n", matrix_size, tile_size, nodes, maxrank, (Ag+uv)/memory_nodes); 
    else
        printf("Matrix size %lld of tile size %d: if maxrank is %d it needs No. of %lf nodes\n", matrix_size, tile_size, maxrank, (Ag+uv)/memory); 

    return 0;
}
