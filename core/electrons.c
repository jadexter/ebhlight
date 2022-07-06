/******************************************************************************
 *                                                                            *
 * ELECTRONS.C                                                                *
 *                                                                            *
 * ELECTRON THERMODYNAMICS                                                    *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"
#include <gsl/gsl_sf_bessel.h>

#if ELECTRONS

void heat_electrons_zone(int i, int j, int k, double Pi[NVAR], double Ps[NVAR],
  double Pf[NVAR], double Dt);
void fixup_electrons_1zone(double P[NVAR]);

void init_electrons()
{
  #if INITELECTRONS
  double uel;

  ZSLOOP(-NG, N1 + NG - 1, -NG, NG + N2 - 1, -NG, NG + N3 - 1) {
    // Set electron internal energy to constant fraction of internal energy
    uel = fel0*P[i][j][k][UU];

    // Initialize entropies
    P[i][j][k][KTOT] = (gam-1.)*P[i][j][k][UU]*pow(P[i][j][k][RHO],-gam);
    P[i][j][k][KEL] = (game-1.)*uel*pow(P[i][j][k][RHO],-game);
  }
  #endif

  bound_prim(P);
}

void heat_electrons(grid_prim_type Pi, grid_prim_type Ps, grid_prim_type Pf,
  double Dt)
{
  timer_start(TIMER_ELECTRON);
  #pragma omp parallel for collapse(3) schedule(dynamic)
  ZLOOP {
    heat_electrons_zone(i, j, k, Pi[i][j][k], Ps[i][j][k], Pf[i][j][k], Dt);
  }
  timer_stop(TIMER_ELECTRON);
}

void heat_electrons_zone(int i, int j, int k, double Pi[NVAR], double Ps[NVAR],
  double Pf[NVAR], double Dt)
{
  double ktotharm, ktotadv, fel;

  // Calculated and advected entropy at final time
  ktotharm = (gam-1.)*Pf[UU]/pow(Pf[RHO],gam);
  ktotadv = Pf[KTOT];

  // Electron heating fraction
  fel = get_fel(i, j, k, Ps);

  // Update final electron entropy according to Ressler+ 2015 Eqn. 27:
  Pf[KEL] += (game-1.)/(gam-1.)*pow(Ps[RHO],gam-game)*fel*(ktotharm-ktotadv);

  // Diagnostics
  struct of_geom *geom = &ggeom[i][j][CENT];
  struct of_state q;
  get_state(Ps, geom, &q);
  double uadv = ktotadv/(gam-1.)*pow(Pf[RHO],gam);
  double Qud = q.ucon[0]*q.ucov[0]*(Pf[UU] - uadv)*pow(Ps[RHO]/Pf[RHO],gam)/Dt;
  // du_e / dtau
  Qvisc_e[i][j][k] = fel*Qud/q.ucov[0];
  // du_p / dtau
  Qvisc_p[i][j][k] = (1-fel)*Qud/q.ucov[0];

  // Reset total entropy
  Pf[KTOT] = ktotharm;
}

// JAD modifying to add alternative e- prescriptions using BETA_HEAT variable
// BETA_HEAT = -gmin from harmpi e.g. 1 is H10, 2 is K18, 3 is W18, 4 is R17
/* adding Kawazura fitting formula equation 2 using felcalc2 as template
   which is supposed to be more accurate than H10 */

 double felhowes(double Trat, double beta)
 {
  double c1, c2, c3, c22, c32, mbeta, qrat, fel;
  mbeta = 2. - 0.2*log10(Trat);
  c1 = 0.92;
  if (Trat <= 1.) {
    c2 = 1.6/Trat;
    c3 = 18. + 5.*log10(Trat);
  } else {
    c2 = 1.2/Trat;
    c3 = 18.;
  }
  c22 = pow(c2, 2.);
  c32 = pow(c3, 2.);
  qrat = c1*(c22+pow(beta,mbeta))/(c32 + pow(beta,mbeta))*exp(-1./beta)*
         pow(MP/ME*Trat,.5);
  fel = 1./(1. + qrat);
  return fel;
}

 double felkawa(double Trat, double beta){
  //    struct of_geom geom ;
  //    double beta, fel, c1, Trat, Tg, Tpr, mbeta, qrat;
  //    double ppres, bsq;
  double qrat, fel;

    qrat = 35./(1.+pow(beta/15.,-1.4)*exp(-0.1/Trat));

    fel = 1./(1.+qrat);

    return fel;
}

/*felwerner:
 ---------
 -- approximate heating fraction from reconnection following Werner+2018 fitting formula
 -- simply varies between fel = 1/4 and 1/2 as a function of sigma_ion
 ***********************************************************************************************/

 double felwerner(double bsq, double rho){
  double fel, sig;
  sig=bsq/rho;
  fel = 1./4.+1./4.*sqrt(sig/5./(2.+sig/5.));
  return fel;
}
/****************************************************************
felrowan():
 ---------
 -- Calculates the fraction of heat given to electrons as a function of beta, sigma following Chael+ based on Rowan+ reconnection simulations
 ***********************************************************************************************/

 double felrowan(double beta, double bsq, double rho, double up, double ue){
  //    struct of_geom geom ;
    double fel, betamax, sigmaw;

    //    get_geometry(i,j,k,CENT,&geom) ;
    //    bsq = bsq_calc(pi[i][j][k],&geom) ; //Magnetic pressure times 2

    /* magnetization parameter */

    //double sig = bsq/pi[i][j][k][RHO];

    //beta = 2.*ppres/(bsq+SMALL);
    //if(beta>1.e20)beta = 1.e20;

    // Chael+2018 equation 9, this version assumes gam=game=gamp OR
    // up >> ue and gam=gamp
    sigmaw = bsq/(rho+gamp*up+game*ue);

    betamax=1./4./sigmaw;
    // Chael+2018 equation 13
    if (beta < betamax) {
      fel=1./2.*exp(-(1.-beta/betamax)/(0.8+sigmaw));
    }
    else {
      fel=1./2.;
    }

    return fel;
}

double get_fel(int i, int j, int k, double P[NVAR])
{
  #if BETA_HEAT == 0
  return fel0;
  #endif

  double fel;
  # if BETA_HEAT > 0
  struct of_geom geom = *get_geometry(i, j, k, CENT);
  double Tpr, uel, Tel, Trat, pres, bsq, beta;

  Tpr = (gam-1.)*P[UU]/P[RHO];
  uel = 1./(game-1.)*P[KEL]*pow(P[RHO],game);
  Tel = (game-1.)*uel/P[RHO];
  if (Tel <= 0.) Tel = SMALL;
  if (Tpr <= 0.) Tpr = SMALL;
  Trat = fabs(Tpr/Tel);
  pres = P[RHO]*Tpr; // Proton pressure
  bsq = bsq_calc(P, &geom);
  beta = pres/bsq*2.;
  if (beta > 1.e20) beta = 1.e20;

  #if BETA_HEAT == 1
  fel = felhowes(Trat,beta);
  #endif

  #if BETA_HEAT==2
  fel = felkawa(Trat,beta);
  #endif

  #if BETA_HEAT==3
  fel = felwerner(bsq,P[RHO]);
  #endif

  #if BETA_HEAT==4
  fel = felrowan(beta,bsq,P[RHO],P[UU],uel);
  #endif

  // enforce a floor on fel for anistropic viscosity
  fel = MY_MAX(fel,0.1);

  #endif

  return fel;

}

// Modified Bessel function of second kind with safe inputs
double safe_Kn(int n, double x)
{
  if (x > 100.) {
    return exp(-x)*sqrt(M_PI/(2.*x));
  } else {
    return gsl_sf_bessel_Kn(n, x);
  }
}

#if RADIATION || COOLING
void coulomb(grid_prim_type Pi, grid_prim_type Ps, grid_prim_type Pf, double Dt)
{
  timer_start(TIMER_ELECTRON);

  #pragma omp parallel for collapse(3)
  ZLOOP {
    double rho = Ps[i][j][k][RHO];
    double thetae = MP/ME*Ps[i][j][k][KEL]*pow(rho,game-1.);
    double ue = Ps[i][j][k][KEL]*pow(rho,game)/(game-1.);
    double up = Ps[i][j][k][UU] - ue;
    double n = rho*Ne_unit;
    double Ti = up*U_unit*(gamp - 1.)/(n*KBOL);
    double thetai = KBOL*Ti/(MP*CL*CL);
    double thetam = 1./(1./thetae + 1./thetai);
    double logCoul = COULOMB_LOG;
    double Te = thetae*ME*CL*CL/KBOL;
    struct of_geom *geom;
    struct of_state q;

    // Sanity checks, although electron fixup routine should catch these
    if (!isnan(Te) && !isnan(Ti) && Te > 0. && Ti > 0.)
    {
      double Qc, term1, term2;

      // Get Coulomb heating rate.
      // Need to handle cases where Thetai < 1e-2, Thetae < 1e-2, and both
      // Thetae and Thetai < 1e-2 separately due to Bessel functions exploding
      double prefac = 3./2.*ME/MP*n*n*logCoul*CL*KBOL*THOMSON*(Ti - Te);
      double thetaCrit = 1.e-2;
      if (thetae < thetaCrit && thetai < thetaCrit) {
        term1 = sqrt(thetam/(M_PI*thetae*thetai/2.));
        term2 = sqrt(thetam/(M_PI*thetae*thetai/2.));
      } else if (thetae < thetaCrit) {
        term1 = exp(-1./thetai)/safe_Kn(2, 1./thetai)*sqrt(thetam/thetae);
        term2 = exp(-1./thetai)/safe_Kn(2, 1./thetai)*sqrt(thetam/thetae);
      } else if (thetai < thetaCrit) {
        term1 = exp(-1./thetae)/safe_Kn(2, 1./thetae)*sqrt(thetam/thetai);
        term2 = exp(-1./thetae)/safe_Kn(2, 1./thetae)*sqrt(thetam/thetai);
      } else {
        term1 = safe_Kn(1, 1./thetam)/(safe_Kn(2, 1./thetae)*safe_Kn(2, 1./thetai));
        term2 = safe_Kn(0, 1./thetam)/(safe_Kn(2, 1./thetae)*safe_Kn(2, 1./thetai));
      }
      term1 *= (2.*pow(thetae + thetai,2) + 1.)/(thetae + thetai);
      term2 *= 2.;
      Qc = prefac*(term1 + term2);

      // Convert to code units
      Qc *= T_unit/U_unit;

      // Update electron internal energy
      geom = &ggeom[i][j][CENT];
      get_state(Ps[i][j][k], geom, &q);

      double ue_f = Pf[i][j][k][KEL]*pow(Pf[i][j][k][RHO],game)/(game-1.);
      ue_f += Qc*Dt/q.ucon[0];

      // Record diagnostic (du_e / dtau)
      Qcoul[i][j][k] = Qc;
      //Qcoul[i][j][k] = q.ucov[0]*Qc;

      // Update electron entropy
      Pf[i][j][k][KEL] = (game-1.)*ue_f*pow(Pf[i][j][k][RHO],-game);
    }
  } // ZLOOP

  timer_stop(TIMER_ELECTRON);
}

#if COOLING
 void electron_cooling(grid_prim_type Ph, double t, double dt)
 {
   // Testing
#pragma omp parallel for collapse(3) schedule(dynamic)
   ZLOOP {
     electron_cooling_zone(i, j, k, Ph[i][j][k], dt);
  }
}

void electron_cooling_zone(int i, int j, int k, double Ph[NVAR], double dt){
  struct of_geom *geom = &ggeom[i][j][CENT];
  struct of_state q;
  double X[NDIM], r, th;
  double thetae, uel, Tel, Y, L, sigma, bsq, Tel_star;
  double Gcov[NDIM];
  // get_state important to update ucov in Gcov!
  get_state(Ph, geom, &q);
  coord(i, j, k, CENT, X);
  bl_coord(X, &r, &th);
  // fprintf(stdout, "coords: %i %i %i %g %g %g %g \n",i,j,k,X[0],X[1],r,th);

  bsq = bsq_calc(Ph, geom);
  sigma = bsq/Ph[RHO];

  double tcool = get_tcool(i, r);

  uel = 1./(game-1.)*Ph[KEL]*pow(Ph[RHO],game);
  // Calculate current electron temperature...
  // Need MP/ME because Ph[KEL] includes *total* (not electron)
  // mass density
  thetae = MP/ME*Ph[KEL]*pow(Ph[RHO],game-1.);
  Tel = thetae*ME*CL*CL/KBOL;
  // Tel = (game-1.)*uel/Ph[RHO]; // JD way...missing MP/ME?

  // calculate cooling rate L following Noble+
  Tel_star = Tel_target*pow(r,-1.*Tel_rslope);
  // fprintf(stdout, "Te target: %.4f\n", Tel_star);
  Y = Tel/Tel_star-1.;
  //L = Omega*uel*pow((Y*(1.+sign(Y))),1./2.);
  //     L = Omega*uel*(Y+abs(Y));
  //     L = Omega*uel*log(1.+(Y+abs(Y)));
  // JD trying very explicit version to avoid any chance of spurious heating from negative L!
  L = 0.;
  //limit Y for now and don't cool sigma > 1
  Y = MY_MIN(Y, 1000.);
  if ((!isnan(Tel)) && (Y > 0.) && (uel > 0.) && (Tel > 0.) && (sigma < 1.)) {
    // this should cool the *electrons* e.g. cooling rate set by their uel not UU
    L = 2.*uel*pow(Y, q_constant)/tcool;
  }
  else {
    L = 0.0;
  }
  // check for supercooling
  if ((uel < dt*L)) {
    fprintf(stderr,"[%i][istart=%i][%i %i %i] supercooling! Y, uel, L: %g %g %g\n",
            mpi_myrank(), global_start[1], i, j, k, Y, uel, L);
    fprintf(stderr, "coords: %g %g %g %g \n",X[0],X[1],r,th);
    // fprintf(stdout, "supercooling! Adjusting tcool to prevent. Y, uel, L: %g %g %g \n", Y, uel, L);
    L = uel/dt - SMALL;
  }
  // AMH added output of L
  Qcool[i][j][k] = L;

  // Implement cooling as a passive sink in local energy conservation (Gcov)
  // update radG, the radiation four-force density (Ryan+ 2015 Eq. 4)
  for (int mu = 0; mu < NDIM; mu++) {
    Gcov[mu] = -L*q.ucov[mu]; // Noble+ 2009 Eqns. 12-13
    radG[i][j][k][mu] = Gcov[mu]*ggeom[i][j][CENT].g;
  }
}

double get_tcool(int i, double r){
  #if TCOOL == 0
  // tcool = constant
  return tcool0;
  #endif

  // Set cooling time according to the orbital time Omega
  // JD: presumably causes problems in any non-BH problem
  #if TCOOL == 1
  double tcool, Omega;
  // If in ISCO, use Noble+ Eq. 16
  if ( r < Risco){
    // Get gcon
    double thetaMidplane = M_PI/2.0;
    double gcon[NDIM][NDIM];
    bl_gcon_func(r, thetaMidplane, gcon);
    // Numerator: g^\phi\mu K_\mu. Contract manually.
    double numerator = gcon[3][0]*Kmu[0]
      + gcon[3][1]*Kmu[1]
      + gcon[3][2]*Kmu[2]
      + gcon[3][3]*Kmu[3];
    // Denominator: g^t\mu K_\mu
    double denominator = gcon[0][0]*Kmu[0]
    + gcon[0][1]*Kmu[1]
    + gcon[0][2]*Kmu[2]
    + gcon[0][3]*Kmu[3];
    Omega = numerator/denominator;
  }
  else{
    // TODO: relativistic correction?
    Omega = 1./(pow(r,3./2.)+a);
  }
  tcool = 1.0/(Omega*tcoolOmega0);
  return tcool;
  #endif

  #if TCOOL == 2
  // Set tcool = A(r-rEH)^p + B
  double p = 1.5; // chosen to imitate 1/Omega dependence
  double B = 1.0; // arbitrarily chosen value of tcool at event horizon
  double Omega0 = 1.0/(pow(Risco, 3./2.) + a); // Omega at the ISCO
  double A = (1.0/Omega0 - B)/(pow(Risco - Reh, p));
  double tcool = A*pow(r - Reh, p) + B;
  // if ( r > Risco){
  return tcool;
  #endif
  #if TCOOL == 3
  // Set tcool = A(r^p-rEH^p) + B
  double p = 1.5; // chosen to imitate 1/Omega dependence
  double B = 1.0; // arbitrarily chosen value of tcool at event horizon
  double Omega0 = 1.0/(pow(Risco, 3./2.) + a); // Omega at the ISCO
  double A = (1.0/Omega0 - B)/(pow(Risco, p) - pow(Reh, p));
  double tcool = A*(pow(r, p) - pow(Reh, p)) + B;
  // if ( r > Risco){
  return tcool;
#endif
}
#endif // COOLING

void apply_rad_force_e(grid_prim_type Prh, grid_prim_type Pr,
  grid_fourvector_type radG, double Dt)
{
  // Apply only to active zones for this proc -- ghost zone four-force
  // depositions already communicated over MPI
  #pragma omp parallel for collapse(3) schedule(dynamic)
  ZLOOP {
    struct of_geom *geom = &ggeom[i][j][CENT];
    struct of_state q;
    double Uel, Urho;
    double C = 0.;

    // Get fluid state at n + 1/2 where radiation four-force is centered
    get_state(Prh[i][j][k], geom, &q);

    for (int mu = 0; mu < NDIM; mu++) {
      C += -q.ucon[mu]*radG[i][j][k][mu];
    }

    // Get fluid state at n+1 for full update
    get_state(Pr[i][j][k], geom, &q);

    // Remove \sqrt{-g} from radG
    C = C/geom->g;

    Urho = Pr[i][j][k][RHO]*q.ucon[0];
    Uel = Pr[i][j][k][KEL]*Urho;

    Uel += Dt*(C*(game-1.)*pow(Prh[i][j][k][RHO],1.-game));

    // Supercooling diagnostics
    if (Uel < 0.) {
      double U_1[NVAR], prim_2[NVAR], U_2[NVAR];
      struct of_state q_1, q_2;
      Nsuper[i][j][k]++;

      // (1) Record total energy density after cooling
      get_state(Pr[i][j][k], &ggeom[i][j][CENT], &q_1);
      primtoflux(Pr[i][j][k], &q_1, 0, &ggeom[i][j][CENT], U_1);

      // (2) Calculate total energy density with zero electron energy density
      PLOOP prim_2[ip] = psupersave[i][j][k][ip];
      double ue = prim_2[KEL]*pow(prim_2[RHO],game)/(game-1.);
      prim_2[UU] -= ue;
      get_state(prim_2, &ggeom[i][j][CENT], &q_2);
      primtoflux(prim_2, &q_2, 0, &ggeom[i][j][CENT], U_2);

      // Subtract (2) from (1); integrated over volume, this is fake energy
      Esuper[i][j][k] += fabs((U_1[UU] - U_2[UU])*dx[1]*dx[2]*dx[3]);
    } // Uel < 0

    Pr[i][j][k][KEL] = Uel/Urho;

    // Reset total entropy
    Pr[i][j][k][KTOT] = (gam-1.)*Pr[i][j][k][UU]*pow(Pr[i][j][k][RHO],-gam);
  } // ZSLOOP
}
#endif // RADIATION

void fixup_electrons(grid_prim_type P)
{
  timer_start(TIMER_ELECTRON);
  #pragma omp parallel for collapse(3) schedule(dynamic)
  ZLOOP {
    fixup_electrons_1zone(P[i][j][k]);
  }
  timer_stop(TIMER_ELECTRON);
}

void fixup_electrons_1zone(double P[NVAR])
{
  double kelmax = P[KTOT]*pow(P[RHO],gam-game)/(tptemin*(gam-1.)/(gamp-1.) + (gam-1.)/(game-1.));
  double kelmin = P[KTOT]*pow(P[RHO],gam-game)/(tptemax*(gam-1.)/(gamp-1.) + (gam-1.)/(game-1.));

  // Replace NANs with cold electrons
  if (isnan(P[KEL])) P[KEL] = kelmin;

  // Enforce maximum Tp/Te
  P[KEL] = MY_MAX(P[KEL], kelmin);

  // Enforce minimum Tp/Te
  P[KEL] = MY_MIN(P[KEL], kelmax);
}
#endif // ELECTRONS
