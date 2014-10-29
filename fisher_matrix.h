#ifndef FISHER_MATRIX3_H
#define FISHER_MATRIX3_H 1

#include <cmath>

const int n_cosmo_param= 4; // omega_m omega_b h ns

const double omega_m= 0.273;
const double omega_b= 0.0456;
const double h= 0.705;
const double sigma_8= 0.812;
const double ns= 0.960;
const double f=std::pow(omega_m, 0.55);

struct Power {
  double k, k_integrated, dk, Pk, dP[n_cosmo_param];
};

struct Volume {
  double dvol, N1, N2, n, n_v, r, z;
};


#endif
