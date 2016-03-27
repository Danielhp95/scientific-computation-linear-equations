#include "gauss.h"


// Returns if Matrix is singular or not.
int reduceToUpperTriangular(double **A, double *y, int N);

//IT IS STUPID BUT START AT A[1][1].
double *Gauss(double **A, double *y, int N) {
  operations = 0;

  int isSingular = reduceToUpperTriangular(A, y, N);
  if (isSingular) {
    return NULL;
  }

  double *roots = malloc(sizeof(double) * (N + 1));

  #pragma omp parallel
  for (int i = N; i >= 0; i--) {
    double known_coefficients = 0.0;
    for (int j = i + 1; j <= N; j++) {
      known_coefficients += roots[j] * A[i][j]; operations += 2;
    }
    roots[i] = (y[i] - known_coefficients) / A[i][i]; operations += 2;
  }
  return roots;
}

int reduceToUpperTriangular(double **A, double *y, int N) {
  #pragma omp parallel 
  for (int i = 1; i < N; i++) {
    for (int j = i + 1; j <= N; j++) {
      if (A[j][i] == 0) {
        // Do nothing
      } else {
        double coefficient = A[j][i] / A[i][i]; operations++;
        for (int t = i; t <= N; t++) {
          A[j][t] -= A[i][t] * coefficient; operations += 2;
        }
        y[j] -= y[i] * coefficient; operations += 2;
        if (A[i][i] == 0.0) {
          return 1; // Matrix has been found to be singular.
        }
      }
    }
  }
  return 0; // Matrix has been succesfully reduced
}

void printMatrix(double **A, int N) {
  for (int i = 1; i <= N; i++) {
    for (int j = 1; j <= N; j++) {
      printf("%.02f ", A[i][j]);
    }
    printf("\n");
  }
}
