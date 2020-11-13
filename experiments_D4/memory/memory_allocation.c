#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>

/*
 * 512 nodes experiments
 * minimal tile size: 1000, and corresponding minrank before 15, maxrank after: 464
 * maximal tile size: 4500, and corresponding minrank before 13, maxrank after: 1079 
 */
int main() {

	struct timeval tstart, tend;
	double t = 0.0;

	double *p;
	int loop = 1000;
#if 0
	int NB = 1000;
	int maxrank = 464; 
	int minrank = 15;
	//int NB = 500;
	//int maxrank = 13; 
	//int minrank = 1;
#else
	int NB = 4500;
	int maxrank = 1079; 
	int minrank = 13;
#endif
	int interval = ( maxrank - minrank ) / 11;

	for( int rank = minrank; rank <= maxrank; rank+=interval ) {
		double flops = 34.0 * NB * rank * rank + 157.0 * rank * rank * rank; 
		double time_dense = flops / ( 1.18 * 1.0e12 / 32 * 0.9 );

		for( int i = 0; i < loop; i++) {
			gettimeofday(&tstart, NULL);
			p = (double *)calloc( 2 * NB * rank, sizeof(double) ); 
			//for( int j = 0; j <  NB * rank; j++ )
			//		p[j] = 1.0;
			free( p );
			gettimeofday(&tend, NULL);
			t += (tend.tv_sec - tstart.tv_sec) + (tend.tv_usec - tstart.tv_usec) / 1000000.0;
		}

		printf("%d %d %d %d %.17le %.17le %.6le\n", NB, maxrank, minrank, rank, t / loop, time_dense, t/loop/time_dense);

	}

	return 0;
}
