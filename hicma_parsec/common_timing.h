#ifndef TIMING_H
#define TIMING_H

#include "parsec/runtime.h"
#include <stdio.h>
#include <sys/time.h>

extern double time_elapsed;
extern double sync_time_elapsed;

#ifdef PARSEC_HAVE_MPI
#define Wtime() MPI_Wtime()
#else
static inline double Wtime(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC_COARSE, &ts);
    return ((double)ts.tv_sec) + (((double)ts.tv_nsec) * 1.0e-9);
}
#endif

#if defined(PARSEC_PROF_TRACE)
#define PARSEC_PROFILING_START() parsec_profiling_start()
#else
#define PARSEC_PROFILING_START()
#endif  /* defined(PARSEC_PROF_TRACE) */

#define TIME_START() do { time_elapsed = Wtime(); } while(0)
#define TIME_STOP() do { time_elapsed = Wtime() - time_elapsed; } while(0)
#define TIME_PRINT(rank, print) do { \
  TIME_STOP(); \
  printf("[%4d] TIME(s) %12.5f : ", rank, time_elapsed); \
  printf print; \
} while(0)

#ifdef PARSEC_HAVE_MPI
# define SYNC_TIME_START() do {                 \
        MPI_Barrier(MPI_COMM_WORLD);            \
        PARSEC_PROFILING_START();                \
        sync_time_elapsed = Wtime();            \
    } while(0)
# define SYNC_TIME_STOP() do {                                  \
        MPI_Barrier(MPI_COMM_WORLD);                            \
        sync_time_elapsed = Wtime() - sync_time_elapsed;        \
    } while(0)
# define SYNC_TIME_PRINT(rank, print) do {                          \
        SYNC_TIME_STOP();                                           \
        if(0 == rank) {                                             \
            printf("[****] TIME(s) %12.5f : ", sync_time_elapsed);        \
            printf print;                                           \
        }                                                           \
  } while(0)

/* overload exit in MPI mode */
#   define exit(ret) MPI_Abort(MPI_COMM_WORLD, ret)

#elif defined(PARSEC_HAVE_LCI)
extern _Noreturn void lci_abort(int exit_code);
# define SYNC_TIME_START() do {                 \
        lc_barrier(*lci_global_ep);             \
        PARSEC_PROFILING_START();               \
        sync_time_elapsed = Wtime();            \
    } while(0)
# define SYNC_TIME_STOP() do {                                  \
        lc_barrier(*lci_global_ep);                             \
        sync_time_elapsed = Wtime() - sync_time_elapsed;        \
    } while(0)
# define SYNC_TIME_PRINT(rank, print) do {                          \
        SYNC_TIME_STOP();                                           \
        if(0 == rank) {                                             \
            printf("[****] TIME(s) %12.5f : ", sync_time_elapsed);        \
            printf print;                                           \
        }                                                           \
  } while(0)

/* overload exit in LCI mode */
#   define exit(ret) lci_abort(ret)

#else
# define SYNC_TIME_START() do { sync_time_elapsed = Wtime(); } while(0)
# define SYNC_TIME_STOP() do { sync_time_elapsed = Wtime() - sync_time_elapsed; } while(0)
# define SYNC_TIME_PRINT(rank, print) do {                           \
        SYNC_TIME_STOP();                                           \
        if(0 == rank) {                                             \
            printf("[****] TIME(s) %12.5f : ", sync_time_elapsed);      \
            printf print;                                           \
        }                                                           \
    } while(0)
#endif

#endif /* TIMING_H */
