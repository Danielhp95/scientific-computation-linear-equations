#include <math.h>
#include "stdio.h"
#include <stdlib.h>
#include <float.h>

double *Gauss(double **, double *, int );
double *BGauss(double **, double *, int , int);
void printMatrix(double **A, int N);
void BprintMatrix(double **A, int N, int B);
