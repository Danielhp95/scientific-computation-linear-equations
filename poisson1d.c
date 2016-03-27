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

double *createYvector(double *y, int N, int value, double delta) {
  y = malloc(sizeof(double) * (N+1));
  for (int i = 1; i <= N; i++) {
    if (i > (N/4) && i <= (N/2)) {
      y[i] = -(delta * delta) * value;
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
    double max_rough, max_smooth, gauss_time, bgauss_time;
    long double gauss_flops, bgauss_flop;
    int rough_x_max, smooth_x_max;
    double *y;
    double **A;
    double delta = 1.0 / (double) N;

    y = createYvector(y, N, ROUGH, delta);
    A = createMatrix(A, N-1);
    A = populateMatrix1D(A, N-1);

    /* ROUGH estimations */
    double start = clock();
    double *rough_roots_gauss = Gauss(A, y, N-1); // Run gauss
    gauss_time = (clock() - start) / CLOCKS_PER_SEC;
    max_rough = max(rough_roots_gauss, N, &rough_x_max);
    gauss_flops = operations * GIGA;

    free(y);
    freeMatrix(A, N-1);
    free(rough_roots_gauss);

    y = createYvector(y, N, ROUGH, delta);
    A = createBandedMatrix(A, N-1);
    start = clock();
    double *roots_bgauss = BGauss(A, y, N-1, 1);  //B = 1 always Run bgauss
    bgauss_time = (clock() - start) / CLOCKS_PER_SEC;
    bgauss_flop = operations * GIGA;

    free(y);
    freeMatrix(A, N-1);
    free(roots_bgauss);

    /* SMOOTH estimationa */
    y = createYvector(y, N, SMOOTH, delta);
    A = createBandedMatrix(A, N-1);
    double *smooth_roots = BGauss(A, y, N-1, 1);
    max_smooth = max(smooth_roots, N, &smooth_x_max);
    free(y);
    freeMatrix(A, N-1);
    free(smooth_roots);

    printf("%i %0.4f %0.4f %0.4f %0.4f %0.8f %0.5Lf %0.8f %0.5Lf\n", N, max_rough, rough_x_max/(double)N, max_smooth, smooth_x_max/(double)N,
                  gauss_time, gauss_flops, bgauss_time, bgauss_flop);
}
