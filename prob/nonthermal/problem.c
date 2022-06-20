/******************************************************************************
 *                                                                            *
 * PROBLEM.C                                                                  *
 *                                                                            *
 * INITIAL CONDITIONS FOR NONTHERMAL TESTING                                  *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"
#include "io.h"

static double plaw;
void set_problem_params()
{
  set_param("plaw", &plaw);
}
void save_problem_params()
{
  WRITE_HDR(plaw, TYPE_DBL);
}

void init_prob()
{
  // double rhoL, PL, u1L, u2L, u3L, B1L, B2L, B3L;
  // double rhoR, PR, u1R, u2R, u3R, B1R, B2R, B3R;


  // tf = 2.;
  
  // rhoL = 1.;
  // PL = 10.;
  // u1L = 1.;
  // u2L = 0.;
  // u3L = 0.;
  // B1L = 0.;
  // B2L = 0.;
  // B3L = 0.;
  // rhoR = 1.;
  // PR = 10.;
  // u1R = 1.0;
  // u2R = 0.;
  // u3R = 0.;
  // B1R = 0.;
  // B2R = 0.;
  // B3R = 0.;
  
  double X[NDIM];

  ZLOOP {
    PLOOP P[i][j][k][ip] = 0.;
    coord(i, j, k, CENT, X);

    // P[i][j][k][RHO] = (X[1] < 0.) ? rhoL : rhoR;
    // P[i][j][k][UU]  = ((X[1] < 0.) ? PL : PR)/(gam - 1.);
    // P[i][j][k][U1]  = (X[1] < 0.) ? u1L : u1R;
    // P[i][j][k][U2]  = (X[1] < 0.) ? u2L : u2R;
    // P[i][j][k][U3]  = (X[1] < 0.) ? u3L : u3R;
    // P[i][j][k][B1]  = (X[1] < 0.) ? B1L : B1R;
    // P[i][j][k][B2]  = (X[1] < 0.) ? B2L : B2R;
    // P[i][j][k][B3]  = (X[1] < 0.) ? B3L : B3R;

    // add some initial gas conditions here if desired

    // inject_nonthermal(P[i][j][k], plaw);
  } // ZLOOP
  // inject_nonthermal(P[NG][NG][NG], plaw);
}

