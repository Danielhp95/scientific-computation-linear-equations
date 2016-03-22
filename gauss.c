#include "gauss.h"


// Returns if Matrix is singular or not.
int reduceToUpperTriangular(double **A, int N);

//IT IS STUPID BUT START AT A[1][1].
double *Gauss(double **A, double *y, int N) {
  // Reduce to upperTriangular


  int isSingular = reduceToUpperTriangular(A, N);
  printf("\nAfter triangulation\n");
  printMatrix(A, N);
  double *roots = malloc(sizeof(double *) * (N+1));
  if (isSingular) {
    return roots;
  } else { // Matrix has a solution
    return roots;
  }
}

int reduceToUpperTriangular(double **A, int N) {
  for (int i = 1; i < N; i++) {
    printf("i: %i \n",i);
    for (int j = i + 1; j <= N; j++) {
      if (A[j][i] == 0) {
        // Do nothing
      } else {
        //A[j][i] = 0;
        double coefficient = A[j][i] / A[i][i];
        for (int t = i; t <= N; t++) {
          printf("t: %i \n", t);
          A[j][t] -= A[i][t] * coefficient;
          printMatrix(A, N);
        }
      }
    }
      //printMatrix(A, N);
  }
}

int isRowSingular(double *row, int N) {
  for (int i = 1; i <=N; i++) {
    if (row[i] != 0) {
      return 0;
    } else {
      break;
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
