//
// Sigma = Cov(delta, u); 2x2 symmetric matrix
//


//
// parameters 0:beta, 1:f*sigma8, 2:rg, 3:sigma_g, 4:sigma_u, 
//            5:omega_m, 6:omega_b, 7:h, 8:ns

class SigmaRedshiftSpace {
 public:
  SigmaRedshiftSpace(const double b_, const double rg_,
		     const double sigma_g_, const double sigma_u_) :
  b(b_), rg(rg_), sigma_g(sigma_g_), sigma_u(sigma_u_), beta(f/b_), beinv(b_/f) {
  }
  int nparam() const { return n_cosmo_param+5; } 
  double get_beta() const { return beta; }
  void operator()(const double k, const double mu, 
		  const double Pdd, const double Pdt, const double Ptt,
		  const double dPdd[], const double dPdt[], const double dPtt[],
		  const double N1, const double N2,
		  double& Psgg, double& Psuu, 
		  double& s11, double& s12, double& s22,
 	          double ds11[], double ds12[], double ds22[]) const {

    //velocity power spectrum in units (100 km/s)^2 (h^{-1} Mpc)^3, or H0=1

    const double kmu= k*mu;
    const double kmu2= kmu*kmu;
    const double mu2= mu*mu;
    const double sigma_g2= sigma_g*sigma_g;
    const double Dg2= 1.0/(1.0 + kmu2*sigma_g2);
    
    const double sinks= sin(k*sigma_u);
    const double cosks= cos(k*sigma_u);
    const double Du= sinks/(k*sigma_u);
    const double Du2= Du*Du;
    const double A= (beinv*beinv + 2.0*rg*beinv*mu2 + mu2*mu2)*f*f*Dg2;
    const double B= mu2/(k*k)*f*f*Du2;
    const double C= mu/k*(rg*beinv + mu2)*f*f*sqrt(Dg2*Du2);

    Psgg= A*Pdd;
    Psuu= B*Ptt;
    s11= A*Pdd + N1;
    s12= C*Pdt;
    s22= B*Ptt + N2;
    
    // 0: beta
    ds11[0]= -2.0*beinv*beinv*(beinv+mu2*rg)*f*f*Dg2*Pdd;
    ds12[0]= -rg*beinv*beinv/(rg*beinv+mu2)*s12;
    ds22[0]= 0.0;

    // 1: f*sigma8
    const double fsigma8= f*sigma_8;
    ds11[1]= 2.0/(fsigma8)*Psgg;
    ds12[1]= 2.0/(fsigma8)*s12;
    ds22[1]= 2.0/(fsigma8)*Psuu;

    // 2: rg
    ds11[2]= 2.0*beinv*mu2*f*f*Dg2*Pdd;
    ds12[2]= beinv/(rg*beinv+mu2)*s12;
    ds22[2]= 0.0;

    // 3: sigma_g
    ds11[3]= -2.0*kmu2*sigma_g*Dg2*Psgg;
    ds12[3]= -kmu2*sigma_g*Dg2*s12;
    ds22[3]= 0.0;

    // 4: sigma_u
    ds11[4]= 0.0;
    ds12[4]= (k*cosks/sinks - 1.0/sigma_u)*s12;
    ds22[4]= 2.0*(k*cosks/sinks - 1.0/sigma_u)*Psuu;

    // 5-8 cosmological parameters
    for(int j=0; j<n_cosmo_param; ++j) {
      ds11[5+j]= A*dPdd[j];
      ds12[5+j]= C*dPdt[j];
      ds22[5+j]= B*dPtt[j];
    }
  }
 private:
  const double b, rg, sigma_g, sigma_u, beta, beinv;
};
