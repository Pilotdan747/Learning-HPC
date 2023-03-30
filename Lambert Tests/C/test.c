#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define DIM1 1000
#define DIM2 1001

int main() {
    double **VinfE = (double**) malloc(sizeof(double)*DIM1); // Gets here
    double **VinfM = (double**) malloc(sizeof(double)*DIM1); // Does not get here -> Seg Fault

    for (int i = 0; i < DIM2; i++) {
        VinfE[i] = (double*) malloc(sizeof(double)*DIM2);
    }

    for (int i = 0; i < DIM1; i++) {
        for (int j = 0; j < DIM2; j++) {
            VinfE[i][j] = i*j + i - j*i;
        }
    }

    double test = fabs(sqrt(-1.0));
    
    printf("Test: %f\n", test);
    printf("Logic: %d\n", test > 0.0001);

    return 0;
}