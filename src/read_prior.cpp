//
// Reads prior from file
// Prior is a 4x4 covariance matrix for Omega_b*h^2, Omega_cdm*h^2, h, and ns
// The covariance matrix is inverted and added to the Fisher matrix


#include "read_prior.h"

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cassert>

#include <gsl/gsl_math.h>
#include <gsl/gsl_linalg.h>

using namespace std;

void read_prior(const char filename[], const int nparam, double* f_prior)
{
  const int n_prior_param= 4;
  const int offset= nparam-n_prior_param;

  FILE* fp= fopen(filename, "r");
  if(fp == 0) {
    cerr << "Unable to read prior file: " << filename << endl;
    throw "read_prior error";
  }

  char buf[128];
  double var[4];
  int i=0;

  gsl_matrix* ccov= gsl_matrix_calloc(n_prior_param, n_prior_param);
  gsl_matrix* ccinv= gsl_matrix_calloc(n_prior_param, n_prior_param);
  gsl_permutation* cperm= gsl_permutation_alloc(n_prior_param);


  while(fgets(buf, 127, fp)) {
    if(buf[0] != '#') {
      int ret= sscanf(buf, "%le %le %le %le ", var, var+1, var+2, var+3);
      assert(ret == 4);

      for(int j=0; j<4; ++j)
	ccov->data[i*n_prior_param+j] = var[j];
      
      i++;
      if(i == 4) break;      
    }
  }

  int signum;
  gsl_linalg_LU_decomp(ccov, cperm, &signum);
  gsl_linalg_LU_invert(ccov, cperm, ccinv);

  for(int i=0; i<n_prior_param; ++i) {
    for(int j=0; j<n_prior_param; ++j) {
      f_prior[(i+offset)*nparam+(j+offset)]= ccinv->data[i*n_prior_param+j];
    }
  }

  fclose(fp);
}
