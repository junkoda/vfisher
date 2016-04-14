#ifndef PRINT_FISHER_MATRIX_H
#define PRINT_FISHER_MATRIX_H 1

#include <gsl/gsl_math.h>
#include <gsl/gsl_linalg.h>

class Printer {
 public:
  virtual ~Printer();
  virtual void print_header(const int nadditional, char const * const names[]);
  virtual void print_k(const double Fij[], const int dim, 
		       const int nadditional, const double additional[]);

};

class Printer1 : public Printer {
 public:
  Printer1(const int nparam, char const * const param_name[], double const param_mean[], int const * const sub_param_index[]);
  virtual ~Printer1();
  virtual void print_header(const int nadditional, char const * const names[]);
  virtual void print_k(const double Fij[], const int dim,
		       const int nadditional, const double additional[]);
 private:
  const int nparam;
  gsl_matrix **cov, **cinv;
  gsl_permutation** perm;

  //const double beta;
  char const * const * const param_name;
  double const * const param_mean;
  int const * const * const sub_param_index;

  void allocate_matrix(const int i);
  void print_sub_matrix_error(const double Fij[], const int nsub, const int isub[]);
};

 


#endif
