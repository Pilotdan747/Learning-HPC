#include <stdio.h>
#include <omp.h>

static long num_steps = 100000000;
double step;

int main() {
    int i; double x, pi, sum = 0.0;

    double t1 = omp_get_wtime();

    step = 1.0/(double) num_steps;

#pragma omp parallel for private(x) reduction(+: sum)
    for (i = 0; i < num_steps; i++) {
        x = (i+0.5)*step;
        sum = sum + 4.0/(1.0 + x*x);
    }

    pi = step*sum;

    double t2 = omp_get_wtime();

    printf("Pi is: %f\n", pi);
    printf("It took %f seconds\n", t2 - t1);
}