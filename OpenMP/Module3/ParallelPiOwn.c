#include <stdio.h>
#include <omp.h>

static long num_steps = 4096/2;
double step;

int main() {
    int j; double x, pi, sum = 0.0;

    double t1 = omp_get_wtime();
    omp_set_num_threads(num_steps);
    printf("Num Threads is: %d\n", omp_get_max_threads());

    step = 1.0/(double) num_steps;

    #pragma omp parallel
    {
        int i = omp_get_thread_num();

        double x = (i+0.5)*step;
        sum = sum + 4.0/(1.0 + x*x);
    }
    pi = step*sum;

    double t2 = omp_get_wtime();

    printf("Pi is: %f\n", pi);
    printf("It took %f seconds\n", t2 - t1);
}