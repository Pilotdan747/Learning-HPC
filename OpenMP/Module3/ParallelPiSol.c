#include <stdio.h>
#include <omp.h>

static long num_steps = 100000000; double step;

// Needed for "False Sharing" protection
#define PAD 16000 // 120 is to PAD out size of Cache Line (128 KB)
#define NUM_THREADS 10 // M1 Pro has 10 cores

int main() {
    int i, nthreads; double pi, sum[NUM_THREADS][PAD];

    double t1 = omp_get_wtime();
    step = 1.0/(double) num_steps;
    omp_set_num_threads(NUM_THREADS);

    #pragma omp parallel
    {
        int i, id, nthrds; double x;

        id = omp_get_thread_num();
        nthrds = omp_get_num_threads();
        
        if (id == 0) nthreads = nthrds;

        for (i = id, sum[id][0] = 0.0; i < num_steps; i = i + nthrds) {
            x = (i+0.5)*step;
            sum[id][0] += 4.0/(1.0 + x*x);
        }
    }

    for (i = 0, pi = 0.0; i < nthreads; i++) pi += sum[i][0]*step;

    double t2 = omp_get_wtime();

    printf("Pi is: %f\n", pi);
    printf("It took %f seconds\n", t2 - t1);
}