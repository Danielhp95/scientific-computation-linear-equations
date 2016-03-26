#include "gauss.h"


// Returns if Matrix is singular or not.
int BreduceToUpperTriangular(double **A, double *y, int N, int B);
int isRowSingular(double *row, int N);

//IT IS STUPID BUT START AT A[1][1].
double *BGauss(double **A, double *y, int N, int B) {

  int isSingular = BreduceToUpperTriangular(A, y, N, B);

  double *roots;
  //roots = (double *) malloc(sizeof(double) * (N + 1));

  int middleIndex = B+1;
  for (int i = N; i <= 1; i--) {
    double known_coefficients = 0.0;
    for (int j = middleIndex; j <= (2*B + 1); j++) {
      if (A[i][j] != 0.0) {
        known_coefficients += roots[i+j] * A[i][j];
      }
      roots[i] = (y[i] - known_coefficients) / A[i][middleIndex];
    }
  }

  return roots;
}

int BreduceToUpperTriangular(double **A, double *y, int N, int B) {
  int middleIndex = B + 1;
  for (int i = 1; i < N; i++) {
    if (A[i][middleIndex] == 0.0) {
      return 1; // Matrix has been found to be singular
    }
    for (int j = 1; (j <= B && j + i <= N); j++) {
      double coefficient = A[i+j][middleIndex-j] / A[i][middleIndex];
      for (int t = 0; t <= B; t++) {
        A[i+j][(middleIndex-j) + t] -= coefficient * A[i][middleIndex + t];
      }
      y[i+j] -= coefficient * y[i];
    }
  }
  return 0; // Matrix has been succesfully reduced
}

int isRowSingular(double *row, int N) {
  for (int i = 1; i <=N; i++) {
    if (row[i] != 0) {
      return 0;
    }
  }
  return 1;
}

void BprintMatrix(double **A, int N, int B) {
  for (int i = 1; i <= N; i++) {
    for (int j = 1; j <= (2*B + 1); j++) {
      printf("%.02f ", A[i][j]);
    }
    printf("\n");
  }
}
