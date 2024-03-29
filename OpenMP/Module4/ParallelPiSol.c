#include <stdio.h>
#include <omp.h>

static long num_steps = 100000000; double step;
#define NUM_THREADS 10

int main() {
    double pi = 0.0;

    double t1 = omp_get_wtime();
    step = 1.0/(double) num_steps;
    omp_set_num_threads(NUM_THREADS);

    #pragma omp parallel
    {
        int i, id, nthrds; double x, sum;

        id = omp_get_thread_num();
        nthrds = omp_get_num_threads();
        
        for (i = id, sum = 0.0; i < num_steps; i = i + nthrds) {
            x = (i + 0.5)*step;
            sum += 4.0/(1.0 + x*x);
        }

        sum = sum*step;

        #pragma omp atomic
            pi += sum;
    }

    double t2 = omp_get_wtime();

    printf("Pi is: %f\n", pi);
    printf("It took %f seconds\n", t2 - t1);
}