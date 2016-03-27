#include "gauss.h"
#include "time.h"

#define ROUGH  80
#define SMOOTH 40
#define GIGA 10e-9

void runApproximation1D(int N);

double max(double *A, int N, int *index) {
  double max = A[1];
  for (int i = 2; i <= N; i++) {
    if (A[i] > max) {
      max = A[i];
      *index = i;
    }
  }
  return max;
}

double **createMatrix(double **A, int N) {
  A = malloc(sizeof(double *) * (N+1));
  for (int i = 0; i <= N; i++) {
    A[i] = malloc(sizeof(double) * (N+1));
    for (int j = 0; j <= N; j++) {
      A[i][j] = 0;
    }
  }
  return A;
}

double **populateMatrix1D(double **A, int N) {

  A[1][1] = -2;
  A[1][2] = 1;

  for (int i = 2; i < N; i++) {
    A[i][i-1] = 1;
    A[i][i] = -2;
    A[i][i+1] = 1;
  }

  A[N][N-1] = 1;
  A[N][N] = -2;
  return A;
}

double **createBandedMatrix(double **A, int N) {
  A = malloc(sizeof(double *) * (N+1));
  for (int i = 1; i <= N; i++) {
    A[i] = malloc(sizeof(double) * (4)); // 2*B + 1 starting at 1 and B = 1.
    A[i][1] = 1;
    A[i][2] = -2;
    A[i][3] = 1;
  }
  A[1][1] = 0;
  A[N][3] = 0;
  return A;
}

double *createYvector(double *y, int N, int value) {
  y = malloc(sizeof(double) * (N+1));
  for (int i = 1; i <= N; i++) {
    if (i > (N/4) && i <= (N/2)) {
      y[i] = value;
    } else {
      y[i] = 0;
    }
  }
  return y;
}

void freeMatrix(double **A, int N) {
  for (int i = 1; i <=N; i++) {
    free(A[i]);
  }
  free(A);
}

int main() {

  /*
      rho(x) is y vector
      poisson_i is roots
       is matrix
       poisson_{i-1} - 2*poisson_i + poisson_{i+1} = -delta^2 * rho_i
       for 0 < i < N
       N = 8, 16, 32,..., 2^30
      (rho(0.25), rho(0.5) = 80 (rough) and 40 (smooth)

       use clock(); for timing
  */

  for (int N = 8; N < (int) pow(2, 30); N *= 2) {
    runApproximation1D(N);
  }


}

void runApproximation1D(int N) {
    double max_rough, max_smooth, gauss_time, gauss_flops,
           bgauss_time, bgauss_flop;
    int rough_x_max, smooth_x_max;
    double *y;
    double **A;

    y = createYvector(y, N, ROUGH);
    A = createMatrix(A, N);
    A = populateMatrix1D(A, N);

    /* ROUGH estimations */
    double start = clock();
    double *rough_roots_gauss = Gauss(A, y, N); // Run gauss
    gauss_time = clock() - start;
    max_rough = max(rough_roots_gauss, N, &rough_x_max);
    gauss_flops = operations * GIGA;

    free(y);
    freeMatrix(A, N);

    y = createYvector(y, N, ROUGH);
    A = createBandedMatrix(A, N);
    start = clock();
    double *roots_bgauss = BGauss(A, y, N, 1);  //B = 1 always Run bgauss
    bgauss_time = clock() - start;
    bgauss_flop = operations * GIGA;

    free(y);
    freeMatrix(A, N);

    /* SMOOTH estimationa */
    y = createYvector(y, N, SMOOTH);
    A = createBandedMatrix(A, N);
    double *smooth_roots = BGauss(A, y, N, 1);
    max_smooth = max(smooth_roots, N, &smooth_x_max);
    free(y);
    freeMatrix(A, N);

    printf("%i %0.4f %i %0.4f %i %0.5f %0.5f %0.5f %0.5f\n", N, max_rough, rough_x_max/N, max_smooth, smooth_x_max/N,
                  gauss_time, gauss_flops, bgauss_time, bgauss_flop);
}

