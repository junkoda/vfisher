//
// Prints marginalised Fisher matrix constraints from given Fisher matrix
//

#include <cstdio>
#include <cassert>
#include <iostream>
#include "fisher_matrix.h"
#include "print_fisher_matrix.h"

using namespace std;

Printer::~Printer()
{

}

void Printer::print_header(const int nadditional, char const * const names[])
{

}

void Printer::print_k(const double Fij[], const int dim, 
		      const int nadditional, const double additional[])
{

}

Printer1::Printer1(const int nparam_, char const * const param_name_[], 
		   const double param_mean_[], 
		   int const * const sub_param_index_[]) :
  nparam(nparam_), param_name(param_name_), param_mean(param_mean_), 
  sub_param_index(sub_param_index_)
{
  // nparam: number of parameters in Fisher Matrix
  // param_name: List parameter names (0 for the end of list)
  // param_mean: value of the parameter (fiducial value)
  // sub_param_index: List of parameter subsets that will be free
  //     (-1 for end of sub indecies, and 0 for end of sub sets)

  cov= (gsl_matrix**) calloc(nparam+1, sizeof(gsl_matrix**));
  cinv= (gsl_matrix**) calloc(nparam+1, sizeof(gsl_matrix**));
  perm= (gsl_permutation**) calloc(nparam+1, sizeof(gsl_permutation**));
}

Printer1::~Printer1()
{
  for(int j=0; j<=nparam; ++j) {
    if(cov[j]) { 
      gsl_permutation_free(perm[j]);
      gsl_matrix_free(cov[j]);
      gsl_matrix_free(cinv[j]);
    }
  }
  
  free(cov);
  free(cinv);
  free(perm);
}

void Printer1::allocate_matrix(const int i)
{
  assert(cov[i] == 0); assert(cinv[i] == 0); assert(perm[i] == 0);
  
  cov[i]= gsl_matrix_calloc(i, i);
  cinv[i]= gsl_matrix_calloc(i, i);
  perm[i]= gsl_permutation_alloc(i);

  assert(cov[i]); assert(cinv[i]); assert(perm[i]);
}

void Printer1::print_sub_matrix_error(const double Fij[], const int nsub, const int isub[])
{
  // Print standard eviation of parameters, marginalized over nsub parameters
  // selected in isub

  assert(nsub <= nparam);
  if(cov[nsub] == 0) {
    allocate_matrix(nsub);
  }

  for(int i=0; i<nsub; ++i) {
    for(int j=0; j<=i; ++j)
      cov[nsub]->data[j*nsub+i] = Fij[isub[j]*nparam+isub[i]];
    for(int j=i+1; j<nsub; ++j)
      cov[nsub]->data[j*nsub+i] = Fij[isub[i]*nparam+isub[j]];
  }

  int signum;
  gsl_linalg_LU_decomp(cov[nsub], perm[nsub], &signum);
  gsl_linalg_LU_invert(cov[nsub], perm[nsub], cinv[nsub]);

  for(int j=0; j<nsub; ++j)
    printf(" %le", sqrt(cinv[nsub]->data[j+j*nsub])/param_mean[isub[j]]);  
}

void Printer1::print_header(const int nadditional, char const * const names_additional[])
{
  for(int j=0; j<nadditional; ++j)
    printf("# %2d %s\n", j+1, names_additional[j]);

  // for ith subset of free parameters
  int icolumn= nadditional+1;
  for(int i=0; sub_param_index[i] != 0; i++) {
    int const * const sub_params= sub_param_index[i];

    // print which parameters are marginalized. Other parameters are fixed.
    printf("# --free: ");
    // for jth parameter in the subset
    for(int j=0; sub_params[j]>=0; ++j)
      printf("%s ", param_name[sub_params[j]]);
    printf("\n");

    for(int j=0; sub_params[j]>=0; ++j)
      printf("# %2d: d%s/%s\n", icolumn++,
	     param_name[sub_params[j]], param_name[sub_params[j]]);
  }
}

//
// Print Fisher Matrix constraints for each k
//
void Printer1::print_k(const double Fij[], const int dim, 
		       const int nadditional, const double additional[])
{
  // parameters: 0:beta, 1:fsigma8, 2:rg, 3:sigma_g, 4:sigma_u,
  //             5:omega_m, 6:omega_b, 7:h

  for(int j=0; j<nadditional; ++j)
    printf("%le ", additional[j]);

  for(int i=0; sub_param_index[i] != 0; i++) {
    int const * const sub_params= sub_param_index[i];
    int nsub= 0;

    for(nsub=0; sub_params[nsub]>=0; ++nsub)
      ;

    assert(nsub > 0);


    print_sub_matrix_error(Fij, nsub, sub_params);
  }

  printf("\n");
}
