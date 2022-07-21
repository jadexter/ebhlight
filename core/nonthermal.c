/******************************************************************************
 *                                                                            *
 * NONTHERMAL.C                                                               *
 *                                                                            *
 * Quantities and functions for nonthermal electron evolution                 *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"

#if NONTHERMAL
void nonthermal_adiab(grid_prim_type Pi, grid_prim_type Pf, double Dt){
    // I worry since the expasion depends on the primitive left and right values that the values may be overwritten as the expansion is calculated...
    // #ifdef ART_ADIAB
    // #define SEMIART_ADIAB (ART_ADIAB)
    // #endif
    #ifdef SEMIART_ADIAB
        double X[NDIM];
        ZLOOPALL{
            coord(i, j, k, CENT, X);
            Pf[i][j][k][U1] = X[1]*SEMIART_ADIAB/3;
            Pf[i][j][k][U2] = X[2]*SEMIART_ADIAB/3;
            Pf[i][j][k][U3] = X[3]*SEMIART_ADIAB/3;
        }
    #endif 
    #pragma omp parallel for collapse(3) schedule(dynamic)
    ZLOOP {
      nonthermal_adiab_zone(i, j, k, Pi, Pf, Dt);
    }
}

void cool_nonthermal(grid_prim_type Pi, grid_prim_type Pf, double Dt){
    #pragma omp parallel for collapse(3) schedule(dynamic)
    ZLOOP {
      cool_nonthermal_zone(Pi[i][j][k], Pf[i][j][k], &(ggeom[i][j][CENT]), Dt);
    }
}

/**
 * @brief Applies radiative cooling to one zone
 * 
 * @param Pr Primitives to calculate the cooling from
 * @param Pf Primitives to apply the cooled result
 * @param geom Geometry in the zone being cooled
 */
void cool_nonthermal_zone(double *Pi, double *Pf, struct of_geom *geom, double Dt){
    double gdot[NTEBINS], ngammas[NTEBINS];
    double dg = log(nteGammas[1])-log(nteGammas[0]);

    NTEGAMMALOOP{
        gdot[ig] = 0;
        ngammas[ig] = Pf[ig+NTESTART];
    } 

    calc_gdot_rad(Pf, geom, gdot);

    // Upwind derivative updates each n(gamma) bin
    NTEGAMMALOOP{
        if(ig == NTEBINS-1){
            Pf[ig+NTESTART] -= Dt*(-gdot[ig]*ngammas[ig])/(dg*nteGammas[ig]);
        }
        else{
            Pf[ig+NTESTART] -= Dt*(gdot[ig+1]*ngammas[ig+1]-gdot[ig]*ngammas[ig])/(dg*nteGammas[ig]);
        }

        if (Pf[ig+NTESTART] < SMALL){
            Pf[ig+NTESTART] = 0;
        }
    }

    // Update thermal electrons with cooling nonthermal electrons
    double ncool = -gdot[0]*ngammas[0];
    double qcool = ME*(nteGammas[0]-1)*ncool;
    double chempot = calc_potential(Pf);
    struct of_state q;

    get_state(Pf, geom, &q);
    apply_thermal_heating(Pf, q, qcool-chempot*ncool,Dt);
}

/**
 * @brief Applies the cooling/heating caused by adiabatic expansion/compression in one zone
 * 
 * @param i Dimension 1 index
 * @param j Dimension 2 index
 * @param k Dimension 3 index
 * @param Pr Full primitives matrix (need neighboring zones as well)
 */
void nonthermal_adiab_zone(int i, int j, int k, grid_prim_type Pi, grid_prim_type Pf, double Dt){
    // TODO: Definitely needs testing...
    // TODO: Include the electrons returning to the thermal distribution (Chael section 3.2 ii)
    // TODO: Viscous dissipation rate and injection terms into thermal and nonthermal pops (Chael section 3.2 iii and eq. 26/29)
    #ifndef ADIABTIC_SCALING
    #define ADIABTIC_SCALING (0)
    #endif

    double adiab = calc_expansion(i,j,k,Pi,Pf,Dt);
    #ifdef ART_ADIAB
    adiab = ART_ADIAB;
    #endif
    double nprime[NTEBINS], deltan[NTEBINS], ngamma[NTEBINS];

    NTEGAMMALOOP ngamma[ig] = Pf[i][j][k][ig+NTESTART];

    // Find the initial values of n and u to compare to the final values
    // double n_tot_start = gamma_integral(ngamma);

    // Find dn/dgamma for use in Chael eq. 47
    nonthermal_adiab_upwind(adiab, ngamma, nprime);

    // Find dtau
    struct of_state q;
    struct of_geom *geom;
    geom = &ggeom[i][j][CENT];
    get_state(Pf[i][j][k], geom, &q);

    double dtau = Dt/(q.ucon[0]);

    #if ADIABTIC_SCALING
        // Find change in n and the resulting real energy change (Chael eq. 47 and 11)
        double ureal[NTEBINS], utot, uexpected[NTEBINS], utotexpected;
    #endif

    NTEGAMMALOOP {
        // Version with semi-analytic derivative
        // deltan[ig] = dtau*(adiab/3.)*( (1+pow(nteGammas[ig],-2))*ngamma[ig] + (nteGammas[ig]-(1/nteGammas[ig]))*nprime[ig] );
        deltan[ig] = dtau*(adiab/3.)*nprime[ig];

        #if ADIABTIC_SCALING
            uexpected[ig] = -1*dtau*ME*(adiab/3.)*(nteGammas[ig]-(1./nteGammas[ig]))*ngamma[ig];
            ureal[ig] = ME*(nteGammas[ig]-1.)*deltan[ig];
        #endif
    }

    #if ADIABTIC_SCALING
        // Rescale change by expected vs actual energies to preserve energy
        utot = gamma_integral(ureal);
        utotexpected = gamma_integral(uexpected);
    #endif

    NTEGAMMALOOP{
        #if ADIABTIC_SCALING
            if ((fabs(utotexpected) > SMALL) && (fabs(utot) > SMALL)) deltan[ig] *= utotexpected/utot;
        #endif

        // Apply the change to the bins
        Pf[i][j][k][ig+NTESTART] += deltan[ig];

        // Floor bins (ie. no negative n)
        if(Pf[i][j][k][ig+NTESTART] < 0){
            Pf[i][j][k][ig+NTESTART] = 0.;
        }
    }
    
    // Catch energy coming out of the nonthermal electrons

    double ncool, qcool;
    // Expanding
    if(adiab > 0){
        ncool = -deltan[0];
        qcool = ncool*ME*(nteGammas[0]-1);
    }
    // Compressing, small amount can escape out the top
    else{
        ncool = -deltan[NTEBINS-1];
        qcool = ncool*ME*(nteGammas[NTEBINS-1]-1);
    }

    double chempot = calc_potential(Pf[i][j][k]);
    apply_thermal_heating(Pf[i][j][k], q, qcool-chempot*ncool,Dt);

}

/**
 * @brief Computes an integral over gamma for a matrix of size NTEBINS
 * 
 * @param ureal Matrix to integrate
 * @return double Integral result
 */
double gamma_integral(double *ureal){
    // TODO: Would logspace or a simpson's integration be better here?
    // TODO: Is the edge handled properly?

    double utot=0;
    double dg;

    NTEGAMMALOOP{
        if(ig == NTEBINS-1)
            dg = pow(10,log10nteGammas[ig]+log10BinSpace) - nteGammas[ig];
        else 
            dg = nteGammas[ig+1] - nteGammas[ig];
        utot += ureal[ig]*dg;
    }
    return utot;
}

/**
 * @brief Performs the partial derivative with respect to gamma of the nonthermal distribution using an upwind finite differencing method
 * 
 * @param adiab four velocity divergence term (adiabtic expansion/contraction term)
 * @param ngamma nonthermal particle number distribution
 * @param nprime matrix to store the resulting derivative
 */
void nonthermal_adiab_upwind(double adiab, double *ngamma, double *nprime)
{
    // TODO: Should probably do some testing here...
    // TODO: Would doing this in logspace be easier? Also confirm it's always upwind

    int sign;
    double upwind, current, dg;
    double fgam[NTEBINS];

    // dg = (log10(NTGAMMAMAX) - log10(NTGAMMAMIN))/((double)(NTEBINS-1))*log(10);

    NTEGAMMALOOP fgam[ig] = (nteGammas[ig]-1/nteGammas[ig])*ngamma[ig];

    NTEGAMMALOOP{
        current = fgam[ig];

        // Expanding
        if(adiab > 0){
            sign = 1;
            if(ig == (NTEBINS-1)){
                upwind = 0;
                dg = (nteGammas[ig]-nteGammas[ig-1]);
            }
            else{
                upwind = fgam[ig+1];
                dg = (nteGammas[ig+1]-nteGammas[ig]);
            } 
        }
        // Compressing
        else{
            sign = -1;
            if(ig == 0){
                upwind = 0;
                dg = (nteGammas[ig+1]-nteGammas[ig]);
            }
            else{
                upwind = fgam[ig-1];
                dg = (nteGammas[ig]-nteGammas[ig-1]);
            } 
        }
        //Derivative
        nprime[ig] = sign*(upwind-current)/dg;
    }

}

/**
 * @brief Calculate the divergence of the four-velocity (expansion/contraction term). Uses {u^\{alpha}}_{;\alpha} = (sqrt(g)*u^\alpha)_,\alpha / sqrt(g)
 * 
 * @param i grid coordinate 1
 * @param j grid coordinate 2
 * @param k grid coordinate 3
 * @param Pr Active primitives
 * @return double {u^\{alpha}}_{;\alpha}
 */
double calc_expansion(int i, int j, int k, grid_prim_type Pi, grid_prim_type Pf, double Dt){
    struct of_state ql, qc, qr;
    struct of_geom *geoml, *geomc, *geomr;
    double Xl[NDIM], Xc[NDIM], Xr[NDIM];
    double du[NDIM];
    double result = 0;


    // Center:
    geomc = &ggeom[i][j][CENT];
    get_state(Pf[i][j][k], geomc, &qc);
    coord(i,j,k,CENT,Xc);


    // Dimension 0:
    double pucon[NDIM];
    ucon_calc(Pi[i][j][k], geomc, pucon);
    du[0] = ( (geomc->g)*(qc.ucon[0])-(geomc->g)*(pucon[0]) )/Dt;

    du[0]=0; // Test if time dimension matters


    // Dimension 1:
    geoml = &ggeom[i-1][j][CENT];
    get_state(Pf[i-1][j][k], geoml, &ql);
    coord(i-1,j,k,CENT,Xl);

    geomr = &ggeom[i+1][j][CENT];
    get_state(Pf[i+1][j][k], geomr, &qr);
    coord(i+1,j,k,CENT,Xr);

    // Could add the option later but I'll just do the center derivative for now
    // du[1] = ( (geomr->g)*(qr.ucon[1])-(geoml->g)*(ql.ucon[1]) )/(Xr[1]-Xl[1]);
    du[1] = ( (geomr->g)*(qr.ucon[1])-(geomc->g)*(qc.ucon[1]) )/(Xr[1]-Xc[1]);


    // Dimension 2:
    geoml = &ggeom[i][j-1][CENT];
    get_state(Pf[i][j-1][k], geoml, &ql);
    coord(i,j-1,k,CENT,Xl);

    geomr = &ggeom[i][j+1][CENT];
    get_state(Pf[i][j+1][k], geomr, &qr);
    coord(i,j+1,k,CENT,Xr);

    // Could add the option later but I'll just do the center derivative for now
    //du[2] = ( (geomr->g)*(qr.ucon[2])-(geoml->g)*(ql.ucon[2]) )/(Xr[2]-Xl[2]);
    du[2] = ( (geomr->g)*(qr.ucon[2])-(geomc->g)*(qc.ucon[2]) )/(Xr[2]-Xc[2]);


    // Dimension 3:
    geoml = &ggeom[i][j][CENT];
    get_state(Pf[i][j][k-1], geoml, &ql);
    coord(i,j,k-1,CENT,Xl);

    geomr = &ggeom[i][j][CENT];
    get_state(Pf[i][j][k+1], geomr, &qr);
    coord(i,j,k+1,CENT,Xr);

    // Could add the option later but I'll just do the center derivative for now
    //du[3] = ( (geomr->g)*(qr.ucon[3])-(geoml->g)*(ql.ucon[3]) )/(Xr[3]-Xl[3]);
    du[3] = ( (geomr->g)*(qr.ucon[3])-(geomc->g)*(qc.ucon[3]) )/(Xr[3]-Xc[3]);


    // Sum and divide:
    DLOOP1 result += du[mu];

    return result/(geomc->g);
}

/**
 * @brief Populate log10nteGammas and nteGammas with appropriate values based on NTGAMMAMAX and NTGAMMAMIN
 * 
 */
void set_nonthermal_gammas()
{
    // TODO: It's not so obvious what I should do about bins. Ie should the values be centered in the bins and I set seperate edge values or something?
    // TODO: May want to use ln instead of log10
    // TODO: when I add injection, there may be some check here to make sure injection min/max is valid

    log10BinSpace = (log10(NTGAMMAMAX) - log10(NTGAMMAMIN))/((double)(NTEBINS-1));
    log10nteGammas[0] = log10(NTGAMMAMIN);
    nteGammas[0] = NTGAMMAMIN;

    double normterms[NTEBINS] = {0};
    double currgamma;

    for (int ig = 1; ig<NTEBINS; ig++){
        log10nteGammas[ig] = log10nteGammas[ig-1] + log10BinSpace;
        nteGammas[ig] = pow(10,log10nteGammas[ig]);

        currgamma = nteGammas[ig];
        if((currgamma>=gammainjmin)&&(currgamma<=gammainjmin)){
            normterms[ig] = ME*(currgamma-1)*pow(currgamma,-PLAW);
        }
    }

    normterm = gamma_integral(normterms);
}

/**
 * @brief Injection of power law distributed nonthermal electrons (see Chael eq. 29)
 * 
 * @param Pr 
 */
void inject_nonthermal(double *Pr, double normalization, double Dt){
    double gammatemp;

    NTEGAMMALOOP{
        gammatemp = nteGammas[ig];
        if((gammatemp <= gammainjmax) && (gammatemp >= gammainjmin)){
            Pr[ig + NTESTART] += Dt*normalization*pow(gammatemp,-PLAW);
        }
    } 
}

void calc_gdot_rad(double *Pr, struct of_geom *geom, double *gdot)
{
    // TODO: Does Bsq direction matter at all?
    // TODO: nion claculation
    // TODO: Inverse compton cooling is just a big ? not sure where to start...

    // Variable Declarations
    struct of_state q;
    #if SYNCHROTRON
    double Bsq;
    #endif
    #if BREMSSTRAHLUNG
    double nion;
    #endif

    get_state(Pr, geom, &q);
    
    #if SYNCHROTRON
    #ifdef ART_SYNCH
    Bsq = pow(ART_SYNCH,2.);
    #else
    Bsq = calc_bsq_cgs(Pr, geom);
    #endif

    //I'm temporarily hard coding in a constant B field to test

    NTEGAMMALOOP gdot[ig] += (-1.292e-11)*Bsq*pow(nteGammas[ig],2); // Eq. 31 Chael + 17 
    #endif

    #if BREMSSTRAHLUNG
    // This just assumes every ion is Hydrogen
    nion = RHO_unit*Pr[RHO]/(MP+ME); // TODO: how to actually find nion... 

    NTEGAMMALOOP gdot[ig] += (-1.37e-16)*nion*nteGammas[ig]*(log(nteGammas[ig])+0.36); // Eq. 32 Chael + 17 
    #endif

    #if INVERSE_COMPTON
    // This one is a bit trickier... Need some help here on how to calulcate Trad and Ehat. See notes for thoughts
    #endif
}

/**
 * @brief Finds the magnitude of the magnetic field in cgs units
 * 
 * @param Pr Primitives in the desired zone
 * @param geom Geomtry in the desired zone
 * @return double B^2 in cgs
 */
double calc_bsq_cgs(double *Pr, struct of_geom *geom){
    double Bp[NDIM], Vcon[NDIM], Ucon[NDIM], Ucov[NDIM], Bcon[NDIM], Bcov[NDIM];
    double Vfac, VdotV, UdotBp;

    Bp[1] = Pr[B1]*B_unit;
    Bp[2] = Pr[B2]*B_unit;
    Bp[3] = Pr[B3]*B_unit;

    Vcon[1] = Pr[U1];
    Vcon[2] = Pr[U2];
    Vcon[3] = Pr[U3];

    // Get Ucov
    VdotV = 0.;
    for(int l = 1; l < NDIM; l++) {
        for(int m = 1; m < NDIM; m++) {
            VdotV += geom->gcov[l][m]*Vcon[l]*Vcon[m];
        }
    }

    Vfac = sqrt(-1./geom->gcon[0][0]*(1. + fabs(VdotV)));
    Ucon[0] = -Vfac*geom->gcon[0][0];
    for(int l = 1; l < NDIM; l++) 
        Ucon[l] = Vcon[l] - Vfac*geom->gcon[0][l];
    lower(Ucon, geom->gcov, Ucov);


    // Get Bcon, Bcov, and B
    UdotBp = 0.;
    for(int l = 1; l < NDIM; l++)
        UdotBp += Ucov[l]*Bp[l];
    Bcon[0] = UdotBp;
    for(int l = 1; l < NDIM; l++)
        Bcon[l] = (Bp[l] + Ucon[l]*UdotBp)/Ucon[0];
    lower(Bcon, geom->gcov, Bcov);

    return dot(Bcon,Bcov);
}

void heat_electrons_zone_nonthermal(int i, int j, int k, double Pi[NVAR], double Ps[NVAR], double Pf[NVAR], double Dt){
    double ktotharm, ktotadv, fel, felth, felnth;
    #ifdef FELNTH
    felnth = FELNTH;
    #else
    felnth = 0.015; // This was the default value from Chael
    #endif

    struct of_geom *geom = &ggeom[i][j][CENT];
    struct of_state qf;
    get_state(Pf, geom, &qf);

    // Calculated and advected entropy at final time
    ktotharm = (gam-1.)*Pf[UU]/pow(Pf[RHO],gam);
    ktotadv = Pf[KTOT];

    // Electron heating fraction
    #ifdef FEL
    fel = FEL;
    #else
    fel = get_fel(i, j, k, Ps);
    #endif
    felth = fel*(1-felnth);

    // Update thermal electron entropy according to Ressler+ 2015 Eqn. 27:
    Pf[KEL] += (game-1.)/(gam-1.)*pow(Ps[RHO],gam-game)*felth*(ktotharm-ktotadv);
    // Update nonthermal electron entropy according to Chael 2017 Eqn. 30:
    double Qtot = (pow(Ps[RHO],gam-1)/(gam-1)) * (Pf[RHO]*(qf.ucon[0])*(ktotharm-ktotadv)/Dt);
    nonthermal_norm = felnth*fel*Qtot/normterm; // normterm is set in set_nonthermal_gammas and is = m_e*int(gam-1)*gam^-p
    inject_nonthermal(Pf, nonthermal_norm, Dt);

    // Diagnostics
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

void apply_thermal_heating(double *Pr, struct of_state q, double heat, double Dt){
    // Update electron internal energy
    double ue_f = Pr[KEL]*pow(Pr[RHO],game)/(game-1.);
    ue_f += heat*Dt/q.ucon[0];

    // Update electron entropy
    Pr[KEL] = (game-1.)*ue_f*pow(Pr[RHO],-game);
}

double calc_potential(double *Pr){
    //From coulomb
    double rho = Pr[RHO];
    double thetae = MP/ME*Pr[KEL]*pow(rho,game-1.);
    double n = rho*Ne_unit;

    return ME * (1. - 0.6*log1p(2.5*thetae) + thetae*(4-1.5*log(thetae*(thetae+0.4)) + log(n)));
}

#endif