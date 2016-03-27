#include "gauss.h"


// Returns if Matrix is singular or not.
int BreduceToUpperTriangular(double **A, double *y, int N, int B);

//IT IS STUPID BUT START AT A[1][1].
double *BGauss(double **A, double *y, int N, int B) {
  operations = 0;

  int isSingular = BreduceToUpperTriangular(A, y, N, B);
  if (isSingular) {
    return NULL;
  }

  double *roots = malloc(sizeof(double) * (N+1));

  int middleIndex = B+1;
  #pragma omp parallel
  for (int i = N; i >= 1; i--) {
    double known_coefficients = 0.0;
    for (int j = 1; j <= B; j++) {
      if (A[i][middleIndex+j] != 0.0) {
        known_coefficients += roots[i+j] * A[i][j+middleIndex]; operations += 2;
      }
    }
    roots[i] = (y[i] - known_coefficients) / A[i][middleIndex]; operations += 2;
  }

  return roots;
}

int BreduceToUpperTriangular(double **A, double *y, int N, int B) {
  int middleIndex = B + 1;
  #pragma omp parallel
  for (int i = 1; i < N; i++) {
    if (A[i][middleIndex] == 0.0) {
      return 1; // Matrix has been found to be singular
    }
    for (int j = 1; (j <= B && j + i <= N); j++) {
      double coefficient = A[i+j][middleIndex-j] / A[i][middleIndex]; operations++;
      for (int t = 0; t <= B; t++) {
        A[i+j][(middleIndex-j) + t] -= coefficient * A[i][middleIndex + t]; operations += 2;
      }
      y[i+j] -= coefficient * y[i]; operations += 2;
    }
  }
  return 0; // Matrix has been succesfully reduced
}

void BprintMatrix(double **A, int N, int B) {
  for (int i = 1; i <= N; i++) {
    for (int j = 1; j <= (2*B + 1); j++) {
      printf("%.02f ", A[i][j]);
    }
    printf("\n");
  }
}
