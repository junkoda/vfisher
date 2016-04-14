#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cassert>
#include <vector>
#include "fisher_matrix.h"
#include "power_spectrum_camb.h"

using namespace std;

inline double power_law_interpolation(const double x, const double x1, const double x2, const double y1, const double y2) {
  const double log_x1= log(x1);
      double theta= (log(x) - log_x1)/(log(x2) - log_x1);
      return exp((1.0-theta)*log(y1) + theta*log(y2));
}

void read_power_spectrum_txt(vector<Power>* const pk, const char filename[], const double kmax)
{
  FILE* fp= fopen(filename, "r");
  if(fp == 0) {
    cerr << "Unable to open power spectrum file: " << filename << endl;
    throw "power_spectrm_camb_error";
  }

  char buf[128]; int i=0;
  while(fgets(buf, 127, fp)) {
    double k, Pk;
    Power pw;
    int ret= sscanf(buf, "%le %le\n", &k, &Pk); assert(ret == 2);
    //(*pk)[i].k= k;
    //(*pk)[i].Pk= Pk;
    pw.k= k; pw.Pk= Pk; assert(k > 0.0);
    pk->push_back(pw);

    //printf("%e %e\n", k, Pk);

    // power-low interpolation for the last point k=kmax
    if(k > kmax) {
      assert(i>0);
      const double log_k0= log((*pk)[i-1].k);
      double theta= (log(kmax) - log_k0)/(log(k) - log_k0);
      (*pk)[i].Pk= exp((1.0-theta)*log((*pk)[i-1].Pk) + theta*log(Pk));
      (*pk)[i].k= kmax;

      //printf("last %e %e\n", (*pk)[i].Pk, (*pk)[i].k);
      break;
    }
    i++;
  }
  fclose(fp);

  if(pk->back().k < kmax) {
    cerr << "Error: Power spectrum file does not have high enough k: " 
	 << pk->back().k << " < " << kmax << endl;
    throw "kmax error";
  }
  
  // dk for trapezoidal integration
  const int n= pk->size(); assert(n > 0);
  double dkprev=0.0;
  for(int j=0; j<n-1; ++j) {
    double dk= (*pk)[j+1].k - (*pk)[j].k;
    //printf("%d dk %e = %e - %e\n", j, dk, (*pk)[j+1].k, (*pk)[j].k);
    //assert(dk > 0.0);
    (*pk)[j].dk= 0.5*(dkprev + dk);
    (*pk)[j].k_integrated= (*pk)[j].k + 0.5*dk;
    dkprev= dk;
  }
  assert(dkprev > 0.0);
  (*pk)[n-1].dk= 0.5*dkprev;
  (*pk)[n-1].k_integrated= (*pk)[n-1].k;

  //cerr << "nk " << pk->size() << endl;
}

void read_power_spectrum_derivative_txt(vector<Power>* const pk, const int iparam, const char filename_plus[], const char filename_minus[], const double dparam)
{
  FILE* fp1= fopen(filename_plus, "r");
  if(fp1 == 0) {
    cerr << "Unable to open power spectrum plus file: " << filename_plus <<endl;
    throw "power_spectrum_plus";
  }

  FILE* fp2= fopen(filename_minus, "r");
  if(fp2 == 0) {
    cerr << "Unable to open power spectrum minus file: " << filename_plus<<endl;
    throw "power_spectrum_minus";
  }
  
  const int n= pk->size(); assert(n>0);
  int ret;
  double k1, Pk1, k2, Pk2;
  double Pk1_prev, Pk2_prev;
  char buf[128];
  for(int j=0; j<n; ++j) {
    fgets(buf, 127, fp1);  
    ret= sscanf(buf, "%le %le", &k1, &Pk1); assert(ret == 2);
    fgets(buf, 127, fp2);
    ret= sscanf(buf, "%le %le", &k2, &Pk2); assert(ret == 2);
    
    const double k= (*pk)[j].k;
    if(j<n-1) {
      assert(abs(k1 - k) < 1.0e-7); assert(abs(k2 - k) < 1.0e-7);
    }
    else {
      assert(j == n-1); assert(j>0); 
      assert(abs(k1 - k2) < 1.0e-7); 
      Pk1= power_law_interpolation(k, (*pk)[j-1].k, k1, Pk1_prev, Pk1);
      Pk2= power_law_interpolation(k, (*pk)[j-1].k, k2, Pk2_prev, Pk2);
    }
    (*pk)[j].dP[iparam] = (Pk1 - Pk2)/(2.0*dparam);

    Pk1_prev= Pk1; Pk2_prev= Pk2;
  }



  fclose(fp1);
  fclose(fp2);
}

vector<Power>* power_spectrum_camb(const double kmax, const char filebase[])
{
  const int nreserve= 500;
  vector<Power>* pk= new vector<Power>;
  pk->reserve(nreserve);
  char filename[256], filename2[256];

  sprintf(filename, "%s0_matterpower.dat", filebase);
  read_power_spectrum_txt(pk, filename, kmax);

  const int nparam= 4;
  const double param[]= {omega_b*h*h, (omega_m-omega_b)*h*h, h, ns};
  for(int j=0; j<nparam; ++j) {
    sprintf(filename, "%s%da_matterpower.dat", filebase, j+1);
    sprintf(filename2, "%s%db_matterpower.dat", filebase, j+1);
    read_power_spectrum_derivative_txt(pk, j, filename, filename2, 0.01*param[j]);
  }
  return pk;
}
