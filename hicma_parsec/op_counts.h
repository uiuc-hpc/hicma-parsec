/**
 * @copyright (c) 2021 King Abdullah University of Science and Technology (KAUST).
 *                     The Universiy of Tennessee and The Universiy of Tennessee Research Foundation.
 *                     All rights reserved.
 * @version 0.1.0
 * @date 2021-01-24
 **/
#ifndef __OP_COUNTS__
#define __OP_COUNTS__
static unsigned long int op_counts(char op, unsigned long int a, unsigned long int b, unsigned long int c, unsigned long int d){
  unsigned long int res = 0;
  if(op == 'q') {//geqrf  if m >= n
    unsigned long int m = a;
    unsigned long int n = b;
    res = 2*m*n*n - (unsigned long int)(2*n*n*n/3.0f) + 2*m*n + (unsigned long int)(17*n/3.0f);
    //printf("%lu %lu %lu\n", m, n, res);
  } 
  else if(op == 'c') {//potrf  
    unsigned long int n = a;
    res = n*n*n/3 - n*n/2.0 + n/6 ;
    //printf("%lu %lu\n",  n, res);
  }
  else if(op == 't') {//trsm  
    unsigned long int m = a;
    unsigned long int n = b;
    int side = c; //1:left 2:right
    if(side == 1)
        res = n*m*m;
    else if(side == 2)
        res = m*n*n;
    else
        fprintf(stderr, "%s %d: invalid side:%d\n", __FILE__, __LINE__, side);
    //printf("%lu %lu\n",  n, res);
  }
  else if(op == 'm') {//gemm  
    unsigned long int m = a;
    unsigned long int n = b;
    unsigned long int k = c;
    res = m*n*k*2;
    //printf("%lu %lu\n",  n, res);
  }
  else if(op == 's') {//svd  
    unsigned long int n = a;
    res = 22*n*n*n;
    //printf("%lu %lu\n",  n, res);
  }
  else if(op == 'o') {//ormqr  
    unsigned long int m = a;
    unsigned long int n = b;
    unsigned long int k = c;
    int side = d; //1:left 2:right
    if(side == 1)
        res = 4*n*m*k-2*n*k*k+3*n*k;
    else if(side == 2)
        res = 4*n*m*k-2*m*k*k+2*m*k+n*k-k*k/2+k/2;
    else
        fprintf(stderr, "%s %d: invalid side:%d\n", __FILE__, __LINE__, side);
    //printf("%lu %lu\n",  n, res);
  }
  else if (op == 'r') {//trmm
    unsigned long int m = a;
    unsigned long int n = b;
    int side = c; //1:left 2:right
    if(side == 1) //left
        res = m*m*n;
    else if(side == 2) //right
        res = m*n*n;
    else
        fprintf(stderr, "%s %d: invalid side:%d\n", __FILE__, __LINE__, side);
  }
  return res;
}

#endif
