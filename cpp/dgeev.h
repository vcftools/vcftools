/*
 * dgeev.h
 * ($Revision: 1 $)
  This file contains a wrapper around the LAPACK dgeev routine,
  used to calculate eigenvalues and associated eigenvectors
  for a square asymmetric matrix H. Since the matrix is asymmetrix,
  both real and imaginary components of the eigenvalues are returned.
  For real, symmetric matricies, the eigenvalues will be real, and
  hence the complex terms equal to zero.

  There are two function calls defined in this header, of the
  forms

    void dgeev(double **H, int n, double *Er, double *Ei)
    void dgeev(double **H, int n, double *Er, double *Ei, double **Evecs)

    H: the n by n matrix that we are solving.
    n: the order of the square matrix H.
    Er: an n-element array to hold the real parts of the eigenvalues.
    Ei: an n-element array to hold the imaginary parts of the eigenvalues.
    Evecs: an n by n matrix to hold the eigenvectors.
*/

#ifndef DGEEV_H_
#define DGEEV_H_

#include <cmath>

void dgeev(double **H, int n, double *Er, double *Ei);
void dgeev(double **H, int n, double *Er, double *Ei, double **Evecs);

double *dgeev_ctof(double **in, int rows, int cols);
void dgeev_ftoc(double *in, double **out, int rows, int cols);
void dgeev_sort(double *Er, double *Ei, int N);
void dgeev_sort(double *Er, double *Ei, double **Evecs, int N);

 
extern "C" void dgeev_(char *jobvl, char *jobvr, int *n, double *a,
		       int *lda, double *wr, double *wi, double *vl,
		       int *ldvl, double *vr, int *ldvr,
		       double *work, int *lwork, int *info);

#endif
