#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include "volume.h"


using namespace std;

// special convension
// ngal_d=0 for zero noise (insted of zero galaxies). same for ngal_v


// Uniform Noise
//  N1 = 1/n_gal_d, and N2= sigma_v^2/n_gal_v are constants
//
vector<Volume>* uniform_noise(const double vol, const double ngal_d, const double ngal_v, const double sigma_v)
{
  vector<Volume>* v= new vector<Volume>;
  const double ngal_d_inv= ngal_d == 0.0 ? 0 : 1.0/ngal_d;
  const double ngal_v_inv= ngal_v == 0.0 ? 0 : 1.0/ngal_v;

  cerr << "ngal_inv " << ngal_d_inv << " " << ngal_v_inv << endl;

  Volume volume;
  volume.dvol= vol;
  volume.N1= 1.0*ngal_d_inv;
  volume.N2= sigma_v*sigma_v*ngal_v_inv;
  //volume.n= ngal;
  //volume.sigma_v= sigma_v;
  v->push_back(volume);

  return v;
}

// Uniform number density for 2 galaxy populations
vector<Volume>* uniform_noise2(const double vol, const double ngal1, const double ngal2)
{
  vector<Volume>* v= new vector<Volume>;
  const double ngal1_inv= ngal1 == 0.0 ? 0 : 1.0/ngal1;
  const double ngal2_inv= ngal2 == 0.0 ? 0 : 1.0/ngal2;

  Volume volume;
  volume.dvol= vol;
  volume.N1= ngal1_inv;
  volume.N2= ngal2_inv;
  v->push_back(volume);

  return v;
}

//
// Uniform Galaxy
//  N1 = 1/n_gal_d is constant, but
//  N2 = sigma_v^2/n_gal_v is fixed fraction of Hubble flow sigma_v= verr*cz
//

std::vector<Volume>* uniform_galaxy(const double steradian, const double ngal_d, const double ngal_v, const double verr_fraction, const double zmax, const int nbin)
{
  vector<Volume>* v= new vector<Volume>;
  v->reserve(nbin);

  const double ngal_d_inv= ngal_d == 0.0 ? 0 : 1.0/ngal_d;
  const double ngal_v_inv= ngal_v == 0.0 ? 0 : 1.0/ngal_v;

  // vol = 1/3*Omega*r_max^3; Omega = field of view [steradian]
  const double r_max= 3000.0*zmax; // c/H0 = 3000.0 [/h Mpc]
  const double vol= steradian/3.0*r_max*r_max*r_max;

  const double sigma_v_star= 3.0; // 300 km/s

  for(int j=0; j<nbin; ++j) {
    Volume volume;
    double r= r_max*double(j+1)/nbin;

    volume.r= r;
    volume.z= zmax*double(j+1)/nbin;
    volume.dvol= vol*(pow(double(j+1)/nbin, 3.0) - pow(double(j)/nbin, 3.0));
    volume.N1= 1.0*ngal_d_inv;

    // H_0^{-1} d = c z = 3.0e5 (km/s) * z
    // sigma_v = verr_fraction*c*z = 3000*verr_fraction*z [100 km/s]

    
    double sigma_v= 3000.0*verr_fraction*zmax*double(j+0.5)/nbin; 
    
    volume.N2= (sigma_v_star*sigma_v_star + sigma_v*sigma_v)*ngal_v_inv;

    v->push_back(volume);
  }

  return v;
}

//
// Reads galaxy number densities n_g(z), n_v(z) from file
//
std::vector<Volume>* nz_list(const char filename[], const double verr_fraction)
{
  FILE* const fp= fopen(filename, "r");
  if(fp == 0) {
    cerr << "Unable to open nz file: " << filename << endl;
    throw "nz_list error";
  }

  printf("# nz_list %s\n", filename);

  const double sigma_v_star= 3.0; // 300 km/s

  vector<Volume>* v= new vector<Volume>;

  char buf[128];
  Volume volume;
  float z, r, dvol, n1, n2;

  while(fgets(buf, 127, fp)) {
    int ret= sscanf(buf, "%e %e %e %e %e", &z, &r, &dvol, &n1, &n2);


    if(ret != 5) {
      cerr << buf;
      cerr << ret << endl;
      break;
    }

    volume.dvol= dvol;
    volume.N1= 1.0/n1;
    volume.n= n1;
    volume.n_v= n2;
    double sigma_v= verr_fraction*r; // internal unit is 100km/s, i.e. H0=1 
    volume.N2= (sigma_v_star*sigma_v_star + sigma_v*sigma_v)/n2;
    volume.z= z;
    
    v->push_back(volume);
  }

  fclose(fp);

  if(v->size() == 0) {
    cerr << "Unable to read nz list lines\n";
    throw "nz_list size 0";
  }

  return v;
}
