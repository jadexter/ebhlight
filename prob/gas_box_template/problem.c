/******************************************************************************
 *                                                                            *
 * PROBLEM.C                                                                  *
 *                                                                            *
 * INITIAL CONDITIONS FOR OPTICALLY THIN BREMSSTRALUNG COOLING                *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"
#include "io.h"

// Problem-specific variables to set at runtime
static double Tp0, Te0, ne;
void set_problem_params()
{
  set_param("Tp0", &Tp0); // Proton temperature in Kelvin
  set_param("Te0", &Te0); // Electron temperature in Kelvin
  set_param("ne", &ne); // total number density
}
void save_problem_params()
{
  WRITE_HDR(Tp0, TYPE_DBL);
  WRITE_HDR(Te0, TYPE_DBL);
  WRITE_HDR(ne, TYPE_DBL);
}

// Initialize dynamical variables
void init_prob()
{
  double rho0 = (MP + ME)*ne;
  double Pp0 = ne*KBOL*Tp0;
  double Pe0 = ne*KBOL*Te0;
  double up0 = Pp0/(gamp-1.);
  double ue0 = Pe0/(game-1.);
  // double ue0 = ne*KBOL*Te0/(game-1.);
  double ug0 = up0 + ue0;
  // Kappa is defined as Pressure/(mass density)^gamma, where
  // mass density is always the total mass density but the pressure
  // is the species' pressure.
  double kappae0 = Pe0/pow(rho0, game)/CL/CL;

  // Alternate definitions.
  // Inspired by void coulomb's definitions.
  // double thetae = KBOL*thetae/(ME*CL*CL);
  // double kappe0 = ME/MP*thetae/pow(rho0, game - 1.0);

  ZLOOP {
    P[i][j][k][RHO] = rho0; // TOTAL mass density (electron + proton)
    P[i][j][j][UU] = ug0; // TOTAL energy density (electron + proton)
    P[i][j][k][KEL] = kappae0;
    P[i][j][k][U1]  = 0.;
    P[i][j][k][U2]  = 0.;
    P[i][j][k][U3]  = 0.;
    P[i][j][k][B1]  = 0.;
    P[i][j][k][B2]  = 0.;
    P[i][j][k][B3]  = 0.;
    fprintf(stdout, "ne, rho0: %g %g\n", ne, rho0);
    fprintf(stdout, "ug0, kappae0: %g %g\n", ug0, kappae0);
  }
}
