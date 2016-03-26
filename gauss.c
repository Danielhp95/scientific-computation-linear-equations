#include "gauss.h"


// Returns if Matrix is singular or not.
int reduceToUpperTriangular(double **A, double *y, int N);
int isRowSingular(double *row, int N);

//IT IS STUPID BUT START AT A[1][1].
double *Gauss(double **A, double *y, int N) {

  int isSingular = reduceToUpperTriangular(A, y, N);

  double *roots;
  //roots = (double *) malloc(sizeof(double) * (N + 1));

  for (int i = N; i >= 0; i--) {
    double known_coefficients = 0.0;
    for (int j = i + 1; j <= N; j++) {
      known_coefficients += roots[j] * A[i][j];
    }
    roots[i] = (y[i] - known_coefficients) / A[i][i];
  }
  for (int i = 1; i <= N; i++) {
    printf("r: %0.2f ", roots[i]);
  }                                                                                                                                                                                                                    
  return roots;
}

int reduceToUpperTriangular(double **A, double *y, int N) {
  for (int i = 1; i < N; i++) {
    for (int j = i + 1; j <= N; j++) {
      if (A[j][i] == 0) {
        // Do nothing
      } else {
        double coefficient = A[j][i] / A[i][i];
        for (int t = i; t <= N; t++) {
          A[j][t] -= A[i][t] * coefficient;
        }
        y[j] -= y[i] * coefficient;
        if (isRowSingular(A[j],N)) {
          return 1; // Matrix has been found to be singular.
        }
      }
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

void printMatrix(double **A, int N) {
  for (int i = 1; i <= N; i++) {
    for (int j = 1; j <= N; j++) {
      printf("%.02f ", A[i][j]);
    }
    printf("\n");
  }
}
