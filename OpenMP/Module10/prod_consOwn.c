/*
**  PROGRAM: A simple serial producer/consumer program
**
**  One function generates (i.e. produces) an array of random values.  
**  A second functions consumes that array and sums it.
**
**  HISTORY: Written by Tim Mattson, April 2007.
*/
#include <omp.h>
#include <stdlib.h>
#include <stdio.h>

#define N        100000000

/* Some random number constants from numerical recipies */
#define SEED       2531
#define RAND_MULT  1366
#define RAND_ADD   150889
#define RAND_MOD   714025
int randy = SEED;
omp_lock_t *lock;
int fillIter;
int sumIter;

/* function to fill an array with random numbers */
void fill_rand(int length, double *a)
{
   int i; 
   for (i=0;i<length;i++) {
     randy = (RAND_MULT * randy + RAND_ADD) % RAND_MOD;
     *(a+i) = ((double) randy)/((double) RAND_MOD);

    if (i%(N/10) == 0){
      #pragma omp flush
      fillIter = i;
    }
   }   

   fillIter = i;
   #pragma omp flush
}

/* function to sum the elements of an array */
double Sum_array(int length, double *a)
{
   int i, count;  double sum = 0.0;

  count = 0;
  while (count < length) {
    for (i=count;i<fillIter;i++) {
      // Need a way to see how far we can go and then wait

      sum += *(a+i);  
    }
    
    count = i;
  }
   return sum; 
}
  
int main()
{
  double *A, sum, runtime;
  int flag = 0;

  A = (double *)malloc(N*sizeof(double));

  runtime = omp_get_wtime();

//omp_init_lock(lock);

#pragma omp parallel
{
  #pragma omp single
    fillIter = 0;

#pragma omp sections
{
    #pragma omp section
      fill_rand(N, A);        // Producer: fill an array of data

    #pragma omp section
      sum = Sum_array(N, A);  // Consumer: sum the array
  }
}
   
  runtime = omp_get_wtime() - runtime;

  printf(" In %f seconds, The sum is %f \n",runtime,sum);
}
 
