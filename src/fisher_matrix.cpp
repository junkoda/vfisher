//
// Multicomponent Fisher Matrix code
// Forecasts cosmological parameter constraints based on
// galaxy number density contrast delta and galaxy peculiar velocity u,
// both in redshift space.
//
// Jun Koda 2011-2014
//

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cassert>
#include <string>
#include <vector>
#include <iomanip>

#include <boost/program_options.hpp>

#include "fisher_matrix.h"
#include "sigma_redshift_space.h"
#include "power_spectrum_camb.h"
#include "volume.h"
#include "print_fisher_matrix.h"
#include "read_prior.h"

using namespace std;
using namespace boost::program_options;


void print_matrix(const double fij[], const int nparam, char const * const param_names[]);

// Galaxy density - peculiar velocity two field Fisher Matrix
struct FisherMatrix2 {
  void operator()(double fij[], const int nparam, 
		    const double d3kbar,
		    const double Psgg, const double Psuu, 
		    const double s11, const double s12, const double s22,
		    const double ds11[], const double ds12[], const double ds22[]
		    ) const {
    double det= s11*s22-s12*s12;
    double detinv2= 1.0/(det*det);

    for(int i=0; i<nparam; ++i) {
      for(int j=i; j<nparam; ++j) {
	fij[nparam*i+j] += 0.5*d3kbar*detinv2*(
	        ( s22*ds11[i] - s12*ds12[i])*( s22*ds11[j] - s12*ds12[j]) +
	        ( s22*ds12[i] - s12*ds22[i])*(-s12*ds11[j] + s11*ds12[j]) +
	        (-s12*ds11[i] + s11*ds12[i])*( s22*ds12[j] - s12*ds22[j]) +
                (-s12*ds12[i] + s11*ds22[i])*(-s12*ds12[j] + s11*ds22[j]));
      }
    }

  }
};

// Galaxy density only
struct FisherMatrixRSD {
  void operator()(double fij[], const int nparam, 
		    const double d3kbar,
		    const double Psgg, const double Psuu, 
		    const double s11, const double s12, const double s22,
		    const double ds11[], const double ds12[], const double ds22[]
		    ) const {
    for(int i=0; i<nparam; ++i) {
      for(int j=i; j<nparam; ++j) {
	fij[nparam*i+j] += 0.5*d3kbar*ds11[i]*ds11[j]/(s11*s11);
      }
    }
  }
};

// Psgu cross power only
struct FisherMatrixCross {
  void operator()(double fij[], const int nparam, 
		    const double d3kbar,
		    const double Psgg, const double Psuu, 
		    const double s11, const double s12, const double s22,
		    const double ds11[], const double ds12[], const double ds22[]
		    ) const {
    const double var= s12*s12 + s11*s22;
    const double fac= var == 0.0 ? 0.0 : d3kbar/var;
    for(int i=0; i<nparam; ++i) {
      for(int j=i; j<nparam; ++j) {
	// 0.5 befor d3kbar is canceled by 0.5 in var[Pgu]
	fij[nparam*i+j] += fac*ds12[i]*ds12[j];
      }
    }
  }
};

template<class Sigma, class FisherMatrix>
void integrate_fisher_matrix2_k(double fij[],
				vector<Power> const * const pdd,
				vector<Power> const * const pdt,
				vector<Power> const * const ptt,
				vector<Volume> const * const volume,
				Sigma sigma,
				FisherMatrix fisher_matrix,
				Printer* const printer)
{
  // Integrate 2-component Fisher Matrix, printing results along k
  //
  //  F_ij = 1/2 Tr[ Sigma^{-1} dSigma/dtheta_i Sigma^{-1} dSigma/dtheta_j ]
  //  Sigma is a 2x2 covariant matrix of (delta u), one of Re or Im part

  const int nmu= 1001;
  const int nk= pdd->size();

  const int nparam= sigma.nparam();

  double Psgg, Psuu;
  double s11, s12, s22;
  double ds11[nparam], ds12[nparam], ds22[nparam];
  
  double nmode= 0.0; // number of independent k-mode (1/2) V_0 Sd^3k/(2pi)^3
  double F_ddamp= 0.0;
  double F_uuonly= 0.0;

  double nmode_eff_g= 0.0, nmode_eff_u= 0.0;

  const int nadditional= 9;
  char const * const names_additional[]= {"k", "nmode", "nmode_eff_g", "nmode_eff_u", "V", "Veff_g", "Veff_d", "dAdd/Add (dbsigma8/bsigma8 Pdd only)", "dAuu/Auu (dfsigma8/fsigma Puu only)"};
  if(printer)
    printer->print_header(nadditional, names_additional);

  for(int ik= 0; ik<nk; ++ik) {
    const double k= (*pdd)[ik].k;
    const double dk= (*pdd)[ik].dk;
    const double Pdd= (*pdd)[ik].Pk;
    double const * const dPdd= (*pdd)[ik].dP;
    const double Pdt= (*pdt)[ik].Pk;
    double const * const dPdt= (*pdt)[ik].dP;
    const double Ptt= (*ptt)[ik].Pk;
    double const * const dPtt= (*ptt)[ik].dP;

    double V= 0.0;
    double Veff_g= 0.0; // spherically (Sdu) averaged Sd3x P/(P+N)
    double Veff_u= 0.0;

    const double eps= 1.0e-2;

    assert(1.0-eps < (*pdt)[ik].k/k && (*pdt)[ik].k/k < 1.0+eps);

    assert(1.0-eps < (*pdt)[ik].dk/dk && (*pdt)[ik].dk/dk < 1.0+eps);
    assert(1.0-eps < (*ptt)[ik].k/k && (*ptt)[ik].k/k < 1.0+eps);
    assert(1.0-eps < (*ptt)[ik].dk/dk && (*ptt)[ik].dk/dk < 1.0+eps);

    for(vector<Volume>::const_iterator 
	  vol= volume->begin(); vol != volume->end(); ++vol) {
      const double N1= vol->N1;
      const double N2= vol->N2;
      for(int imu=0; imu<nmu; ++imu) {
	const double mu= (double)imu/(nmu-1);
	const double dmu= (1.0-0.5*((imu==0)||(imu==(nmu-1))))/(nmu-1);
	// Trapezoidal integration. dmu=0.5 for first and last step
	//   (f_0 + f_1)/2*dmu + (f_1 + f_2)/2*dmu + ...
	//   = f_0*dmu/2 + f_1*dmu + ...

	const double fac= 1.0/(2.0*M_PI*M_PI)*dk*dmu; 
	// d3kbar= S d^3k/(2pi)^3 = 1/(2*pi)^3 * 2pi * S k^2 dk 2*S_0^1 dmu
	//       = fac*k*k

	const double d3kbar= vol->dvol*fac*k*k;
	nmode += 0.5*d3kbar;

	sigma(k, mu, Pdd, Pdt, Ptt, dPdd, dPdt, dPtt, N1, N2,
	      Psgg, Psuu, s11, s12, s22, ds11, ds12, ds22);

	const double s22_inv= s22 > 0.0 ? 1.0/s22 : 0.0; // Psuu = 0 for mu=0
	const double var_Psuu= (Psuu*s22_inv)*(Psuu*s22_inv);

	F_ddamp += 2.0*d3kbar*(Psgg/s11)*(Psgg/s11);
	F_uuonly += 2.0*d3kbar*var_Psuu;

	V      += vol->dvol*dmu;
	Veff_g += vol->dvol*dmu*(Psgg/s11)*(Psgg/s11);
	Veff_u += vol->dvol*dmu*var_Psuu;
	
	nmode_eff_g += 0.5*d3kbar*(Psgg/s11)*(Psgg/s11);	
	nmode_eff_u += 0.5*d3kbar*var_Psuu;
	
	// S = Sigma = Cov(delta, u)
	// F_ij= (1/2)Tr[ S^{-1} dS/dtheta_i S^{1-} dS/dtheta_j ]
	fisher_matrix(fij, nparam, d3kbar,
		      Psgg, Psuu, s11, s12, s22, ds11, ds12, ds22);
		      


      }
    }
    // print at each k
    double additional[]= {(*pdd)[ik].k_integrated, 
			  nmode, nmode_eff_g, nmode_eff_u, 
			  V, Veff_g, Veff_u, 
			  1.0/sqrt(F_ddamp), 1.0/sqrt(F_uuonly)};

    if(printer)
      printer->print_k(fij, nparam, nadditional, additional);
  }
}

template<class Sigma, class FisherMatrix>
void integrate_fisher_matrix2_z(double* fij,
				vector<Power> const * const pdd,
				vector<Power> const * const pdt,
				vector<Power> const * const ptt,
				vector<Volume> const * const volume,
				Sigma sigma,
				FisherMatrix fisher_matrix,
				Printer* const printer)
{
  // Integrate 2-component Fisher Matrix, printing results along z

  const int nmu= 1001;
  const int nk= pdd->size();

  const int nparam= sigma.nparam();

  double Psgg, Psuu;
  double s11, s12, s22;
  double ds11[nparam], ds12[nparam], ds22[nparam];
  
  double nmode= 0.0; // number of independent k-mode (1/2) SdV Sd^3k/(2pi)^3
  double F_ddamp= 0.0;
  double F_uuonly= 0.0;

  const int nadditional= 7;
  char const * const names_additional[]= {"z", "V", "nmode", "nmode_eff_g", "nmode_eff_u", "dAdd/Add (dbsigma8/bsigma8 Pdd only)", "dAuu/Auu (dfsigma8/fsigma Puu only)"};
  if(printer)
    printer->print_header(nadditional, names_additional);

  double V= 0.0;
  double nmode_eff_g= 0.0;
  double nmode_eff_u= 0.0;

  for(vector<Volume>::const_iterator 
	vol= volume->begin(); vol != volume->end(); ++vol) {

    V += vol->dvol;

    for(int ik= 0; ik<nk; ++ik) {
     const double k= (*pdd)[ik].k;
      const double dk= (*pdd)[ik].dk;
      const double Pdd= (*pdd)[ik].Pk;
      double const * const dPdd= (*pdd)[ik].dP;
      const double Pdt= (*pdt)[ik].Pk;
      double const * const dPdt= (*pdt)[ik].dP;
      const double Ptt= (*ptt)[ik].Pk;
      double const * const dPtt= (*ptt)[ik].dP;
      
      const double eps= 1.0e-2;
      assert(1.0-eps < (*pdt)[ik].k/k && (*pdt)[ik].k/k < 1.0+eps);
      assert(1.0-eps < (*pdt)[ik].dk/dk && (*pdt)[ik].dk/dk < 1.0+eps);
      assert(1.0-eps < (*ptt)[ik].k/k && (*ptt)[ik].k/k < 1.0+eps);
      assert(1.0-eps < (*ptt)[ik].dk/dk && (*ptt)[ik].dk/dk < 1.0+eps);

      const double N1= vol->N1; //1.0/vol->n;
      const double N2= vol->N2; //vol->sigma_v*vol->sigma_v/vol->n;
      for(int imu=0; imu<nmu; ++imu) {
	const double mu= (double)imu/(nmu-1);
	const double dmu= (1.0-0.5*((imu==0)||(imu==(nmu-1))))/(nmu-1);
	// Trapezoidal integration. dmu=0.5 for first and last step
	//   (f_0 + f_1)/2*dmu + (f_1 + f_2)/2*dmu + ...
	//   = f_0*dmu/2 + f_1*dmu + ...

	const double fac= 1.0/(2.0*M_PI*M_PI)*dk*dmu; 
	// d3kbar= S d^3k/(2pi)^3 = 1/(2*pi)^3 * 2pi * S k^2 dk 2*S_0^1 dmu
	//       = fac*k*k
	const double d3kbar= vol->dvol*fac*k*k;
	nmode += 0.5*d3kbar;

	sigma(k, mu, Pdd, Pdt, Ptt, dPdd, dPdt, dPtt, N1, N2,
	      Psgg, Psuu, s11, s12, s22, ds11, ds12, ds22);

	const double s22_inv= s22 > 0.0 ? 1.0/s22 : 0.0;
	const double var_Psuu= (Psuu*s22_inv)*(Psuu*s22_inv);

	F_ddamp += 2.0*d3kbar*(Psgg/s11)*(Psgg/s11);
	F_uuonly += 2.0*d3kbar*var_Psuu;

	nmode_eff_g += 0.5*d3kbar*(Psgg/s11)*(Psgg/s11);	
	nmode_eff_u += 0.5*d3kbar*var_Psuu;

	// S = Sigma = Cov(delta, u)
	// F_ij= (1/2)Tr[ S^{-1} dS/dtheta_i S^{1-} dS/dtheta_j ]
	fisher_matrix(fij, nparam, d3kbar,
		      Psgg, Psuu, s11, s12, s22, ds11, ds12, ds22);	
      }
    }
    // print at each k
    double additional[]= {vol->z, V, nmode, nmode_eff_g, nmode_eff_u, 
			  1.0/sqrt(F_ddamp), 1.0/sqrt(F_uuonly)};

    if(printer)
      printer->print_k(fij, nparam, nadditional, additional);
  }
}

template<class Sigma, class FisherMatrix>
void integrate_fisher_matrix2_dk(double fij[],
				vector<Power> const * const pdd,
				vector<Power> const * const pdt,
				vector<Power> const * const ptt,
				vector<Volume> const * const volume,
				Sigma sigma,
				FisherMatrix fisher_matrix,
				Printer* const printer,
				const double dk_out)
{
  // Integrate 2-component Fisher Matrix, printing results for kbin
  //
  //  F_ij = 1/2 Tr[ Sigma^{-1} dSigma/dtheta_i Sigma^{-1} dSigma/dtheta_j ]
  //  Sigma is a 2x2 covariant matrix of (delta u), one of Re or Im part

  const int nmu= 1001;
  const int nk= pdd->size();

  const int nparam= sigma.nparam();

  double Psgg, Psuu;
  double s11, s12, s22;
  double ds11[nparam], ds12[nparam], ds22[nparam];
  
  double nmode= 0.0; // number of independent k-mode (1/2) V_0 Sd^3k/(2pi)^3
  double F_ddamp= 0.0;
  double F_uuonly= 0.0;
  double F_ddbetalimit= 0.0; // dbeta/beta cosmic variance limit for Pdd only.

  double nmode_eff_g= 0.0, nmode_eff_u= 0.0;

  const int nadditional= 9;
  char const * const names_additional[]= {"k", "nmode", "nmode_eff_g", "nmode_eff_u", "V", "Veff_g", "Veff_d", "dAdd/Add (dbsigma8/bsigma8 Pdd only)", "dAuu/Auu (dfsigma8/fsigma Puu only)"};
  if(printer)
    printer->print_header(nadditional, names_additional);

  // save Fij prior
  double* fij_init= (double*) calloc(sizeof(double), nparam*nparam);
  for(int j=0; j<nparam*nparam; ++j)
    fij_init[j]= fij[j];

  //const double dk_out= 0.01;
  double k_out= dk_out;

  for(int ik= 0; ik<nk; ++ik) {
    const double k= (*pdd)[ik].k;
    double dk= (*pdd)[ik].dk;
    const double Pdd= (*pdd)[ik].Pk;
    double const * const dPdd= (*pdd)[ik].dP;
    const double Pdt= (*pdt)[ik].Pk;
    double const * const dPdt= (*pdt)[ik].dP;
    const double Ptt= (*ptt)[ik].Pk;
    double const * const dPtt= (*ptt)[ik].dP;

    double V= 0.0;
    double Veff_g= 0.0; // spherically (Sdu) averaged Sd3x P/(P+N)
    double Veff_u= 0.0;

    const double eps= 1.0e-2;

    assert(1.0-eps < (*pdt)[ik].k/k && (*pdt)[ik].k/k < 1.0+eps);

    assert(1.0-eps < (*pdt)[ik].dk/dk && (*pdt)[ik].dk/dk < 1.0+eps);
    assert(1.0-eps < (*ptt)[ik].k/k && (*ptt)[ik].k/k < 1.0+eps);
    assert(1.0-eps < (*ptt)[ik].dk/dk && (*ptt)[ik].dk/dk < 1.0+eps);

    for(vector<Volume>::const_iterator 
	  vol= volume->begin(); vol != volume->end(); ++vol) {
      const double N1= vol->N1;
      const double N2= vol->N2;
      for(int imu=0; imu<nmu; ++imu) {
	const double mu= (double)imu/(nmu-1);
	const double dmu= (1.0-0.5*((imu==0)||(imu==(nmu-1))))/(nmu-1);
	// Trapezoidal integration. dmu=0.5 for first and last step
	//   (f_0 + f_1)/2*dmu + (f_1 + f_2)/2*dmu + ...
	//   = f_0*dmu/2 + f_1*dmu + ...

	const double fac= 1.0/(2.0*M_PI*M_PI)*dk*dmu; 
	// d3kbar= S d^3k/(2pi)^3 = 1/(2*pi)^3 * 2pi * S k^2 dk 2*S_0^1 dmu
	//       = fac*k*k

	const double d3kbar= vol->dvol*fac*k*k;
	nmode += 0.5*d3kbar;

	sigma(k, mu, Pdd, Pdt, Ptt, dPdd, dPdt, dPtt, N1, N2,
	      Psgg, Psuu, s11, s12, s22, ds11, ds12, ds22);

	const double s22_inv= s22 > 0.0 ? 1.0/s22 : 0.0; // Psuu = 0 for mu=0
	const double var_Psuu= (Psuu*s22_inv)*(Psuu*s22_inv);

	F_ddamp += 2.0*d3kbar*(Psgg/s11)*(Psgg/s11);
	F_uuonly += 2.0*d3kbar*var_Psuu;

	V      += vol->dvol*dmu;
	Veff_g += vol->dvol*dmu*(Psgg/s11)*(Psgg/s11);
	Veff_u += vol->dvol*dmu*var_Psuu;
	
	nmode_eff_g += 0.5*d3kbar*(Psgg/s11)*(Psgg/s11);	
	nmode_eff_u += 0.5*d3kbar*var_Psuu;
	
	// S = Sigma = Cov(delta, u)
	// F_ij= (1/2)Tr[ S^{-1} dS/dtheta_i S^{1-} dS/dtheta_j ]
	fisher_matrix(fij, nparam, d3kbar,
		      Psgg, Psuu, s11, s12, s22, ds11, ds12, ds22);
		      


      }
    }
    // print at k_out
    if((*pdd)[ik].k_integrated > k_out) {
      double additional[]= {(*pdd)[ik].k_integrated, 
			    nmode, nmode_eff_g, nmode_eff_u, 
			    V, Veff_g, Veff_u, 
			    1.0/sqrt(F_ddamp), 1.0/sqrt(F_uuonly)};
      
      if(printer)
	printer->print_k(fij, nparam, nadditional, additional);

      // reset
      for(int j=0; j<nparam*nparam; ++j)
	fij[j]= fij_init[j];
      nmode= 0.0; F_ddamp= 0.0; F_uuonly= 0.0; F_ddbetalimit= 0.0; 
      nmode_eff_g= 0.0; nmode_eff_u= 0.0;

       
      k_out += dk_out;
    }
  }

  free(fij_init);
}

int main(int argc, char* argv[])
{
  //
  // program_options
  //
  options_description opt("fisher_matrix3 [options] mode");
  opt.add_options()
    ("help,h", "display this help")
    ("mode", "k, z, dk, matrix, final")
    ("noise", value<string>()->default_value("nz"), "uniform, uniform_ngal, or <nz filename>")
    ("ngal", value<double>()->default_value(1.0e-3, "1.0e-3"), "uniform galaxy number density [(1/h Mpc)^-3]")
    ("ngalv", value<double>()->default_value(0.0, "ngal"), "ngal with peculiar velocity")
    ("sigmav", value<double>()->default_value(6000.0, "6000"), "uniform velocity error [km/s] for noise=uniform")
    ("steradian", value<string>()->default_value("4pi"), "field of view")
    ("verr", value<double>()->default_value(0.2, "0.2"), "fractional error in peculiar velocity for noise=uniform_ngal or nz")
    ("pk", value<string>()->default_value("camb"), "power pectra camb, spt, or rpt")
    ("kmax", value<double>()->default_value(0.3, "0.3"), "kmax [h Mpc^-1]")
    ("zmax", value<double>()->default_value(0.1, "0.1"), "zmax")
    ("nrbin", value<int>()->default_value(200), "number of volume bins for noise=unifrom or uniform_ngal")
    ("b", value<double>()->default_value(1.0, "1"), "bias")
    ("rg", value<double>()->default_value(1.0, "1"), "cross correlation factor")
    ("sigma_g", value<double>()->default_value(3.0, "3"), "damping parameter of Pgg")
    ("sigma_u", value<double>()->default_value(13.0, "13"), "damping parameter of Puu")
    ("rsd-only", "redshift space distortion (Pgg) only")
    ("cross-only", "cross power (Pgu) only")
    ("planck-prior", value<string>()->default_value("planck_cov_nolens_8pw.txt"), "Add planck prior Fisher Matrix")
    ("dk-out", value<double>()->default_value(0.02, "0.02"), "dk_out in dk mode")
    ("pk-dir", value<string>()->default_value("power_spectra"), "directory containing power spectra")
    ;

  positional_options_description p;
  p.add("mode", -1);
  variables_map vm;
  store(command_line_parser(argc, argv).
	options(opt).positional(p).run(), vm);
  notify(vm);    
  
  if(vm.count("help") || ! vm.count("mode")) {
    cout.setf(ios::scientific);
    cout << setiosflags(ios::scientific) << opt << "\n"
	 << "mode:\n"
            "  k             print uncertainties along k\n"
            "  z             print uncertainties along z\n"
            "  matrix        print Fisher Matrix\n";

    return 0;
  }

  const string mode= vm["mode"].as<string>();
  const string noise= vm["noise"].as<string>();

  const double kmax= vm["kmax"].as<double>();

  // Survey Parameters
  const double ngal= vm["ngal"].as<double>();
  const double ngalv= vm["ngalv"].as<double>() == 0.0 ? 
                      ngal : vm["ngalv"].as<double>();
  const string str= vm["steradian"].as<string>();
  const double steradian = str.compare(str.size()-2, 2, "pi") == 0 ?
    M_PI*atof(str.substr(0, str.size()-2).c_str()) : atof(str.c_str());

  const double zmax= vm["zmax"].as<double>();
  const double verr_fraction= vm["verr"].as<double>();

  // Numerical Precision Parameters (neglected for --noise=<filename>
  const int nrbin= vm["nrbin"].as<int>();

  // Model Fiducial Values
  const double b= vm["b"].as<double>();
  const double rg= vm["rg"].as<double>();
  const double sigma_g= vm["sigma_g"].as<double>();
  const double sigma_u= vm["sigma_u"].as<double>();

  const string power_spectrum= vm["pk"].as<string>();

  // velocity unit is 100 km/s

  printf("# vfisher\n");
  printf("# kmax zmax %f %f\n", kmax, zmax);

  //
  // Power Spectrum selection
  //
  vector<Power> *pdd= 0, *pdt= 0, *ptt= 0;

  string pk_dir= vm["pk-dir"].as<string>();

  if(power_spectrum == "camb") {
    printf("# power_spectrum camb (linear)\n");
    string filebase= pk_dir + "/linear/camb";

    pdd= power_spectrum_camb(kmax, filebase.c_str());
    pdt= pdd;
    ptt= pdd;
  }
  else if(power_spectrum == "spt") {
    printf("# power_spectrum 3pt\n");
    const string filebase_dd= pk_dir + "/3pt/pkd/camb";
    pdd= power_spectrum_camb(kmax, filebase_dd.c_str());
    
    string filebase_dt= pk_dir + "/3pt/pkdt/camb";
    pdt= power_spectrum_camb(kmax, filebase_dt.c_str());

    string filebase_tt= pk_dir + "/3pt/pkt/camb";
    ptt= power_spectrum_camb(kmax, filebase_tt.c_str());
  }
  else if(power_spectrum == "rpt") {
    printf("# power_spectrum rpt\n");
    string filebase_dd= pk_dir + "/rpt/pdd/camb";
    pdd= power_spectrum_camb(kmax, filebase_dd.c_str());
    
    string filebase_dt= pk_dir + "/rpt/pdt/camb";
    pdt= power_spectrum_camb(kmax, filebase_dt.c_str());

    string filebase_tt= pk_dir + "/rpt/ptt/camb";
    ptt= power_spectrum_camb(kmax, filebase_tt.c_str());
  }
  else {
    cerr << "Unknown power spectrum: " << power_spectrum << endl;
    return 1;
  }
  
  //
  // Volume and Noise term selection
  //
  vector<Volume>* volume= 0;
  if(noise == "uniform") {
    // uniform velocity noise
    const double sigma_v= vm["sigmav"].as<double>() / 100.0; 
                                              // internal velocity unit 100km/s
    const double rmax= 3000.0*zmax;
    const double vol= steradian/3.0*rmax*rmax*rmax;
    volume= uniform_noise(vol, ngal, ngalv, sigma_v);
    printf("# uniform_noise vol ngal ngalv sigma_v %e %e %e %e\n", 
	   vol, ngal, ngalv, sigma_v);
  }
  else if(noise == "uniform_ngal") {
    // velocity noise proportional to distance
    volume= uniform_galaxy(steradian,ngal, ngalv, verr_fraction, zmax, nrbin);
    assert(volume->size() > 0);
    printf("# uniform_galaxy steradian ngal ngalv verr_fraction zmax nrbin %e %e %e %f %f %d\n", 
	   steradian, ngal, ngalv, verr_fraction, zmax, nrbin);
  }
  else {
    // reads galaxy number densitys from file
    volume= nz_list(noise.c_str(), verr_fraction);
    printf("# nz_file %s verr_fraction nrbin %e %d\n", 
	   noise.c_str(), verr_fraction, (int)volume->size());
  }

  // Result printing options
  char const * const param_names[]= {"beta", "fsigma8", "rg", "sigma_g", "sigma_u", "omegabh2", "omegach2", "h", "ns"};
  const double param_values[] = {f/b, f*sigma_8, rg, sigma_g, sigma_u, omega_b*h*h, omega_m*h*h, h, ns};

  const int sub_params_beta[]= {0, -1};
  const int sub_params_fsigma8[]= {1, -1};
  const int sub_params2[]= {0,1, -1};       // free parameters: beta fsigma8
  const int sub_params3[]= {0,1,2, -1};     // free parameters: beta fsigma8 rg
  const int sub_params4[]= {0,1,3,4, -1};   //  beta fsigma8 sigma_g sigma_u
  const int sub_params5[]= {0,1,2,3,4, -1}; //  beta fsigma8 rg sigma_g sigma_u
  const int sub_paramsc[]= {0,1,5,6,7,8, -1}; // beta fsigma8 cosmo params
  const int sub_params9[]= {0,1,2,3,4,5,6,7,8, -1}; // all

  int const * const sub_indeces_2[]= 
    {sub_params_beta, sub_params_fsigma8, sub_params2, sub_params3, sub_params4, sub_params5, sub_paramsc, sub_params9, 0};


  const int sub_params_rsd3[]= {0,1,3, -1};     //  beta fsigma8 sigma_g
  const int sub_params_rsd4[]= {0,1,2,3, -1};   //  beta fsigma8 rg sigma_g
  const int sub_params_rsdc[]= {0,1,5,6,7,8, -1}; // all except sigma_u
  const int sub_params_rsd8[]= {0,1,2,3,5,6,7,8, -1}; // all except sigma_u
  
  int const * const sub_indeces_rsdonly[]=
    {sub_params_beta, sub_params_fsigma8, sub_params2, sub_params3, sub_params_rsd3, sub_params_rsd4, sub_params_rsdc, sub_params_rsd8, 0};



  SigmaRedshiftSpace sigma(b, rg, sigma_g, sigma_u);
  Printer* printer= 0;

  printf("# model sigma_redshift_space b rg sigma_g sigma_u %f %f %f %f\n",
	 b, rg, sigma_g, sigma_u);


  const int nparam= sigma.nparam();
  double* fij= (double*) calloc(sizeof(double), nparam*nparam);

  { // 100% prior on sigma_g and sigma_u
    const int isigma_g= 3;
    fij[nparam*isigma_g + isigma_g] = 1.0/(sigma_g*sigma_g);
    const int isigma_u= 4;
    fij[nparam*isigma_u + isigma_u] = 1.0/(sigma_u*sigma_u);
  }

  // Read PLANCK or any prior from file

  if(vm.count("planck-prior")) {
    const string prior_filename= vm["planck-prior"].as<string>();
    printf("# planck prior file %s\n", prior_filename.c_str());
    read_prior(prior_filename.c_str(), nparam, fij);
  }

  if(vm.count("rsd-only")) {
    printf("# --rsd-only\n");
    // This is a case using redshift sample only "RSD only"
    // No sigma_u in RSD-only

    printer= new Printer1(nparam, param_names, param_values, 
			  sub_indeces_rsdonly);

    if(mode == "k") {
      printf("# mode k\n");
      integrate_fisher_matrix2_k(fij, pdd, pdt, ptt, volume, sigma, 
				 FisherMatrixRSD(), printer);
    }
    else if(mode == "z") {
      printf("# mode z\n");
      integrate_fisher_matrix2_z(fij, pdd, pdt, ptt, volume, sigma, 
				 FisherMatrixRSD(), printer);    
    }
    else if(mode == "matrix") {
      printf("# mode matrix\n");
      integrate_fisher_matrix2_k(fij, pdd, pdt, ptt, volume, sigma, 
				 FisherMatrixRSD(), 0);
      print_matrix(fij, nparam, param_names);
    }
    else if(mode == "final") {
      printf("# mode final\n");
      integrate_fisher_matrix2_k(fij, pdd, pdt, ptt, volume, sigma, 
				 FisherMatrixRSD(), 0);

      if(printer)
	printer->print_k(fij, nparam, 0, 0);
    }
    else if(mode == "dk") {
      const double dk= vm["dk-out"].as<double>(); assert(dk > 0.0);
      printf("# mode dk= %f\n", dk);
      integrate_fisher_matrix2_dk(fij, pdd, pdt, ptt, volume, sigma, 
				  FisherMatrixRSD(), printer, dk);
    }
    else {
      cerr << "Unknown mode: " << mode << endl;
      return 1;
    }
  }
  else if(vm.count("cross-only")) {
    printf("# --cross-only\n");
    printer= 
      new Printer1(nparam, param_names, param_values, sub_indeces_2);

    if(mode == "k") {
      printf("# mode k\n");
      integrate_fisher_matrix2_k(fij, pdd, pdt, ptt, volume, sigma, 
				 FisherMatrixCross(), printer);
    }
    else if(mode == "z") {
      printf("# mode z\n");
      integrate_fisher_matrix2_z(fij, pdd, pdt, ptt, volume, sigma, 
				 FisherMatrixCross(), printer);    
    }
    else if(mode == "matrix") {
      printf("# mode matrix\n");
      integrate_fisher_matrix2_k(fij, pdd, pdt, ptt, volume, sigma, 
				 FisherMatrixCross(), 0);
      print_matrix(fij, nparam, param_names);
    }
    else if(mode == "final") {
      printf("# mode final\n");
      integrate_fisher_matrix2_k(fij, pdd, pdt, ptt, volume, sigma, 
				 FisherMatrixCross(), 0);

      if(printer)
	printer->print_k(fij, nparam, 0, 0);
    }
    else if(mode == "dk") {
      const double dk= vm["dk-out"].as<double>(); assert(dk > 0.0);
      printf("# mode dk= %f\n", dk);
      integrate_fisher_matrix2_dk(fij, pdd, pdt, ptt, volume, sigma, 
				 FisherMatrixCross(), printer, dk);
    }
    else {
      cerr << "Unknown mode: " << mode << endl;
      return 1;
    }
  }
  else {
    printf("# two-field\n");
    printer= 
      new Printer1(nparam, param_names, param_values, sub_indeces_2);
    if(mode == "k") {
      printf("# mode k\n");
      integrate_fisher_matrix2_k(fij, pdd, pdt, ptt, volume, sigma, 
				 FisherMatrix2(), printer);
    }
    else if(mode == "z") {
      printf("# mode z\n");
      integrate_fisher_matrix2_z(fij, pdd, pdt, ptt, volume, sigma, 
				 FisherMatrix2(), printer);    
    }
    else if(mode == "matrix") {
      printf("# mode matrix\n");
      integrate_fisher_matrix2_k(fij, pdd, pdt, ptt, volume, sigma, 
				 FisherMatrix2(), 0);
      print_matrix(fij, nparam, param_names);
    }
    else if(mode == "final") {
      printf("# mode final\n");
      integrate_fisher_matrix2_k(fij, pdd, pdt, ptt, volume, sigma, 
				 FisherMatrix2(), 0);

      if(printer)
	printer->print_k(fij, nparam, 0, 0);
    }
    else if(mode == "dk") {
      const double dk= vm["dk-out"].as<double>(); assert(dk > 0.0);
      printf("# mode dk= %f\n", dk);
      printf("# mode k\n");
      integrate_fisher_matrix2_dk(fij, pdd, pdt, ptt, volume, sigma, 
				 FisherMatrix2(), printer, dk);
    }
    else {
      cerr << "Unknown mode: " << mode << endl;
      return 1;
    }
  }


  return 0;
}

void print_matrix(const double fij[], const int nparam, char const * const param_names[])
{
  printf("# ");
  for(int j=0; j<nparam; ++j)
    printf("%s ", param_names[j]);
  printf("\n");

  for(int i=0; i<nparam; ++i) {
    for(int j=0; j<i; ++j)
      printf("%e ", fij[j*nparam+i]);
    for(int j=i; j<nparam; ++j)
      printf("%e ", fij[i*nparam+j]);

    printf("\n");
  }
}
