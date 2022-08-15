/******************************************************************************
 *                                                                            *
 * NONTHERMAL_SCATTERING.C                                                    *
 *                                                                            *
 * NONTHERMAL ELECTRON MODIFICATION TO COMPTON SCATTERING                     *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"

#if RADIATION
#if NONTHERMAL

#define NW 220
#define NT 90
#define MINW (1.e-16)
#define MAXW (1.e10)
#define MINT (0.001)
#define MAXT (1.e11)
#define HOTCROSS "hotcross.dat"

double nt_alpha_inv_scatt(double nu, const struct of_microphysics *m, double *P);
double nt_compton_cross_calc(double epsilon, double *ngamma);
void nt_dNdgammae(double *ngamma, double *ndot, double *dgamma);
double nt_hc_klein_nishina(double we);
double nt_get_scatt_bias(double nu, const struct of_microphysics *m, double uph, double *P);
void sample_nonthermal(double *ngamma, double *gamma_e, double *beta_e);

double nt_alpha_inv_scatt(double nu, const struct of_microphysics *m, double *P){
  double ngamma[NTEBINS];
  NTEGAMMALOOP ngamma[ig] = P[ig+NTESTART];
  
  double Eg = HPL*nu/(ME*CL*CL);
  double kappa = nt_compton_cross_calc(Eg, ngamma);

  return nu*kappa*(m->Ne);
}

double nt_get_scatt_bias(double nu, const struct of_microphysics *m, double uph, double *P){
    double Thetae = m->Thetae;
    double amp = 1. + 4.*Thetae - 2.*pow(Thetae,3./2.) + 16.*pow(Thetae,2.);
    double bias = tune_scatt*amp;

    // Ensure bias is in reasonable physical bounds
    if (bias < 1.) {
        bias = 1.;
    }

    // Another BOUND_BIAS method. Assumes large hotspots, may work fine instead
    // assuming ~GM/c^2 length scale for hot spots.
    double dl = Rout_rad*L_unit;
    double dtau = ( (alpha_inv_scatt(nu, m)+nt_alpha_inv_scatt(nu, m, P)) /nu)*dl;
    bias = MY_MIN(bias, 1./dtau);
    bias = MY_MAX(bias, 1.);

    return bias;
}


#define NTDMUE (0.05)
 
double nt_compton_cross_calc(double epsilon, double *ngamma)
{
    double cross, gammae, betae, epsilon_e;
    double ndot[NTEBINS], dgammae[NTEBINS];
    nt_dNdgammae(ngamma, ndot, dgammae);

    if(isnan(epsilon)) {
    fprintf(stderr, "NAN Compton cross section: %g\n", epsilon);
    return 0.;
    }

    double dmue = NTDMUE;

    // Integrate over mu_e and gamma_e, where mu_e is the cosine of the angle
    // between K and U_e, and the angle k is assumed to lie, wlog, along the z
    // z axis
    cross = 0.;
    // cross = integral over gamma and mu of the functions given in Dolence eq. 27
    for (double mue = -1. + 0.5*dmue; mue < 1.; mue += dmue) {
        NTEGAMMALOOP {
            gammae = nteGammas[ig];
            betae = sqrt(gammae*gammae-1)/gammae;
            epsilon_e = epsilon*gammae*(1-mue*betae);

            cross += dmue*dgammae[ig]*(1-mue*betae)*nt_hc_klein_nishina(epsilon_e)*ndot[ig];

            if(isnan(cross)) {
                fprintf(stderr, "NAN cross section: %g %g %g %g %g %g\n", epsilon, mue, gammae, ndot[ig], nt_hc_klein_nishina(epsilon_e),(1-mue*betae));
            }
        }
    }

    return cross*THOMSON;
}

void nt_dNdgammae(double *ngamma, double *ndot, double *dgamma){
    NTEGAMMALOOP{
        if(ig != NTEBINS-1){
            dgamma[ig] = (nteGammas[ig+1]-nteGammas[ig]);
            ndot[ig] = (ngamma[ig+1]-ngamma[ig])/dgamma[ig];
        }
        else{
            dgamma[ig] = (pow(10,log10nteGammas[ig]+log10BinSpace)-nteGammas[ig]);
            ndot[ig] = ndot[ig-1]-(ngamma[ig]/dgamma[ig]);
        }
    }
}

double nt_hc_klein_nishina(double we)
{
  double sigma;

  if (we < 1.e-3) return(1. - 2.*we);

  sigma = (3./4.)*(
    2./(we*we) +
    (1./(2.*we) -
    (1. + we)/(we*we*we))*log(1. + 2.*we) +
    (1. + we)/((1. + 2.*we)*(1. + 2.*we))
    );

  return sigma;
}

void sample_nonthermal(double *ngamma, double *gamma_e, double *beta_e)
{
    double weight_sum = 0;
    NTEGAMMALOOP{
        weight_sum += ngamma[ig];
    }

    double rnd = get_rand()*weight_sum;
    int bin = 0;
    NTEGAMMALOOP{
        if(rnd<ngamma[ig]){
            bin = ig;
            break;
        }
        rnd -= ngamma[ig];
    }
    if (bin==NTEBINS-1){
        *gamma_e = nteGammas[bin];
    }
    else{
        *gamma_e = nteGammas[bin]+get_rand()*(nteGammas[bin+1]-nteGammas[bin]);
    }
    *beta_e = sqrt(1. - 1./((*gamma_e)*(*gamma_e)));
}

int nt_sample_electron(double k[NDIM], double p[NDIM], double Thetae, double Ne, double *Prad)
{
  double beta_e, mu, phi, cphi, sphi, gamma_e, sigma_KN;
  double K, sth, cth, x1, n0dotv0, v0, v1;
  double n0x, n0y, n0z;
  double v0x, v0y, v0z;
  double v1x, v1y, v1z;
  double v2x, v2y, v2z;
  int sample_cnt = 0;
  double ngamma[NTEBINS];
  NTEGAMMALOOP ngamma[ig] = Prad[ig+NTESTART];
  double n_enth = gamma_integral(ngamma);
  int sampledNT = 0;

  do {
    sampledNT = 0;
    if (get_rand() < (Ne-n_enth)/Ne){
        sample_beta(Thetae, &gamma_e, &beta_e);
    }
    else{
        sample_nonthermal(ngamma, &gamma_e, &beta_e);
        sampledNT = 1;
    }
    
    mu = sample_mu(beta_e);

    // Sometimes |mu| > 1 from roundoff error. Fix it
    if (mu > 1.) mu = 1.;
    else if (mu < -1.) mu = -1;

    // Frequency in electron rest frame
    K = gamma_e*(1. - beta_e*mu)*k[0];

    // Avoid numerical craziness for small K
    if (K < 1.e-3) {
      sigma_KN = 1. - 2.*K + 5.2*K*K - 13.3*K*K*K + 1144*K*K*K*K/35.;
    } else{
      // Klein-Nishina cross-section / Thomson
      sigma_KN = (3./(4.*K*K))*(2. + K*K*(1. + K)/
        ((1. + 2.*K)*(1. + 2.*K)) + (K*K - 2.*K - 2.)/(2.*K)*log(1. + 2.*K));
    }

    x1 = get_rand();

    sample_cnt++;

    if(sample_cnt > 1000000) {
      fprintf(stderr,"in sample_electron mu, gamma_e, K, sigma_KN, x1: %g %g %g %g %g %g\n",
        Thetae, mu, gamma_e, K, sigma_KN, x1);

      // Kluge to prevent stalling for large values of \Theta_e
      Thetae *= 0.5 ;
      sample_cnt = 0 ;
    }
  } while (x1 >= sigma_KN);

  // First unit vector for coordinate system
  v0x = k[1];
  v0y = k[2];
  v0z = k[3];
  v0 = sqrt(v0x*v0x + v0y*v0y + v0z*v0z);
  v0x /= v0;
  v0y /= v0;
  v0z /= v0;

  // Pick zero-angle for coordinate system
  get_ran_dir_3d(&n0x, &n0y, &n0z);
  n0dotv0 = v0x*n0x + v0y*n0y + v0z*n0z;

  // Second unit vector
  v1x = n0x - (n0dotv0)*v0x;
  v1y = n0y - (n0dotv0)*v0y;
  v1z = n0z - (n0dotv0)*v0z;

  // Normalize
  v1 = sqrt(v1x*v1x + v1y*v1y + v1z*v1z);
  v1x /= v1;
  v1y /= v1;
  v1z /= v1;

  // Find one more unit vector using cross product; automatically normalized
  v2x = v0y*v1z - v0z*v1y;
  v2y = v0z*v1x - v0x*v1z;
  v2z = v0x*v1y - v0y*v1x;

  // Resolve new momentum vector along unit vectors and create a four-vector p
  phi = get_rand()*2.*M_PI; // uniform orientation
  sphi = sin(phi);
  cphi = cos(phi);
  cth = mu;
  sth = sqrt(1. - mu*mu);

  p[0] = gamma_e;
  p[1] = gamma_e*beta_e*(cth*v0x + sth*(cphi*v1x + sphi*v2x));
  p[2] = gamma_e*beta_e*(cth*v0y + sth*(cphi*v1y + sphi*v2y));
  p[3] = gamma_e*beta_e*(cth*v0z + sth*(cphi*v1z + sphi*v2z));

  if (beta_e < 0) {
    fprintf(stderr, "betae error: %g %g %g %g\n", p[0], p[1], p[2], p[3]);
  }

  return sampledNT;
}


#endif // NONTHERMAL
#endif // RADIATION