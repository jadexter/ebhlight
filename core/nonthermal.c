/******************************************************************************
 *                                                                            *
 * RADIATION.C                                                                *
 *                                                                            *
 * MODEL-INDEPENDENT RADIATION QUANTITIES                                     *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"

#define SYNCHROTRON (1)
#define BREMSSTRAHLUNG (0)
#define INVERSE_COMPTON (0)

#if NONTHERMAL
#endif
void step_nonthermal(grid_prim_type Pr)
{
    #pragma omp parallel for collapse(3)
    ZLOOP{
        // inject_nonthermal(Pr[i][j][k],-2.5);
        // nonthermal_adiab(i, j, k, Pr);
        // cool_nonthermal(Pr[i][j][k], &(ggeom[i][j][CENT]));
    }
}

/**
 * @brief Applies the cooling/heating caused by adiabatic expansion/compression in one zone
 * 
 * @param i Dimension 1 index
 * @param j Dimension 2 index
 * @param k Dimension 3 index
 * @param Pr Full primitives matrix (need neighboring zones as well)
 */
void nonthermal_adiab(int i, int j, int k, grid_prim_type Pr){
    // TODO: Definitely needs testing...
    // TODO: Include the electrons returning to the thermal distribution (Chael section 3.2 ii)
    // TODO: Viscous dissipation rate and injection terms into thermal and nonthermal pops (Chael section 3.2 iii and eq. 26/29)

    double adiab = calc_expansion(i,j,k,Pr);
    adiab = -5e-3;
    double nprime[NTEBINS], deltan[NTEBINS], ngamma[NTEBINS];

    NTEGAMMALOOP ngamma[ig] = Pr[i][j][k][ig+NTESTART];

    // Find the initial values of n and u to compare to the final values
    // double n_tot_start = gamma_integral(ngamma);

    // Find dn/dgamma for use in Chael eq. 47
    nonthermal_adiab_upwind(adiab, ngamma, nprime);

    // Find dtau
    struct of_state q;
    struct of_geom *geom;
    geom = &ggeom[i][j][CENT];
    get_state(Pr[i][j][k], geom, &q);

    double dtau = dt/(q.ucon[0]);
    

    #ifdef ADIABTIC_SCALING
        // Find change in n and the resulting real energy change (Chael eq. 47 and 11)
        double ureal[NTEBINS], utot, uexpected[NTEBINS], utotexpected;
    #endif

    NTEGAMMALOOP {
        // Version with semi-analytic derivative
        // deltan[ig] = dtau*(adiab/3.)*( (1+pow(nteGammas[ig],-2))*ngamma[ig] + (nteGammas[ig]-(1/nteGammas[ig]))*nprime[ig] );
        deltan[ig] = dtau*(adiab/3.)*nprime[ig];

        #ifdef ADIABTIC_SCALING
            uexpected[ig] = -1*dtau*ME*(adiab/3.)*(nteGammas[ig]-(1./nteGammas[ig]))*ngamma[ig];
            ureal[ig] = ME*(nteGammas[ig]-1.)*deltan[ig];
        #endif
    }

    #ifdef ADIABTIC_SCALING
        // Rescale change by expected vs actual energies to preserve energy
        utot = gamma_integral(ureal);
        utotexpected = gamma_integral(uexpected);
    #endif

    NTEGAMMALOOP{
        #ifdef ADIABTIC_SCALING
            if ((fabs(utotexpected) > SMALL) && (fabs(utot) > SMALL)) deltan[ig] *= utotexpected/utot;
        #endif

        // Apply the change to the bins
        Pr[i][j][k][ig+NTESTART] += deltan[ig];

        // Floor bins (ie. no negative n)
        if(Pr[i][j][k][ig+NTESTART] < 0){
            Pr[i][j][k][ig+NTESTART] = 0.;
        }
    }
    
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
double calc_expansion(int i, int j, int k, grid_prim_type Pr){
    struct of_state ql, qc, qr;
    struct of_geom *geoml, *geomc, *geomr;
    double Xl[NDIM], Xc[NDIM], Xr[NDIM];
    double du[NDIM];
    double result = 0;


    // Center:
    geomc = &ggeom[i][j][CENT];
    get_state(Pr[i][j][k], geomc, &qc);
    coord(i,j,k,CENT,Xc);


    // Dimension 0:
    du[0]=0;


    // Dimension 1:
    geoml = &ggeom[i-1][j][CENT];
    get_state(Pr[i-1][j][k], geoml, &ql);
    coord(i-1,j,k,CENT,Xl);

    geomr = &ggeom[i+1][j][CENT];
    get_state(Pr[i+1][j][k], geomr, &qr);
    coord(i+1,j,k,CENT,Xr);

    // Could add the option later but I'll just do the center derivative for now
    du[1] = ( (geomr->g)*(qr.ucon[1])-(geoml->g)*(ql.ucon[1]) )/(Xr[1]-Xl[1]);


    // Dimension 2:
    geoml = &ggeom[i][j-1][CENT];
    get_state(Pr[i][j-1][k], geoml, &ql);
    coord(i,j-1,k,CENT,Xl);

    geomr = &ggeom[i][j+1][CENT];
    get_state(Pr[i][j+1][k], geomr, &qr);
    coord(i,j+1,k,CENT,Xr);

    // Could add the option later but I'll just do the center derivative for now
    du[2] = ( (geomr->g)*(qr.ucon[2])-(geoml->g)*(ql.ucon[2]) )/(Xr[2]-Xl[2]);


    // Dimension 3:
    geoml = &ggeom[i][j][CENT];
    get_state(Pr[i][j][k-1], geoml, &ql);
    coord(i,j,k-1,CENT,Xl);

    geomr = &ggeom[i][j][CENT];
    get_state(Pr[i][j][k+1], geomr, &qr);
    coord(i,j,k+1,CENT,Xr);

    // Could add the option later but I'll just do the center derivative for now
    du[3] = ( (geomr->g)*(qr.ucon[3])-(geoml->g)*(ql.ucon[3]) )/(Xr[3]-Xl[3]);


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

    for (int ig = 1; ig<NTEBINS; ig++){
        log10nteGammas[ig] = log10nteGammas[ig-1] + log10BinSpace;
        nteGammas[ig] = pow(10,log10nteGammas[ig]);
    }
}

/**
 * @brief Applies radiative cooling to one zone
 * 
 * @param Pr Primitives for the zone to cool
 * @param geom Geometry in the zone being cooled
 */
void cool_nonthermal(double *Pr, struct of_geom *geom){
    double gdot[NTEBINS], ngammas[NTEBINS];
    double dg = log(nteGammas[1])-log(nteGammas[0]);

    NTEGAMMALOOP{
        gdot[ig] = 0;
        ngammas[ig] = Pr[ig+NTESTART];
    } 

    calc_gdot_rad(Pr, geom, gdot);

    // Upwind derivative updates each n(gamma) bin
    NTEGAMMALOOP{
        if(ig == NTEBINS-1){
            Pr[ig+NTESTART] -= dt*(-gdot[ig]*ngammas[ig])/(dg*nteGammas[ig]);
        }
        else{
            Pr[ig+NTESTART] -= dt*(gdot[ig+1]*ngammas[ig+1]-gdot[ig]*ngammas[ig])/(dg*nteGammas[ig]);
        }

        if (Pr[ig+NTESTART] < SMALL){
            Pr[ig+NTESTART] = 0;
        }
    }
}

/**
 * @brief Injection of power law distributed nonthermal electrons (see Chael eq. 29)
 * 
 * @param Pr 
 */
void inject_nonthermal(double *Pr, double powerlaw){
    // TODO: Add parameters passed at compile time for maximum and minimum injection gammas, C and p
    // TODO: This currently only handles constant, user defined injection but this should definitely change...

    double gammainjmax = 1e5;
    double gammainjmin = 500;
    double normalization = 1;
    double gammatemp;

    NTEGAMMALOOP{
        gammatemp = nteGammas[ig];
        if((gammatemp <= gammainjmax) && (gammatemp >= gammainjmin)){
            Pr[ig + NTESTART] += dt*normalization*pow(gammatemp,powerlaw);
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
    Bsq = calc_bsq_cgs(Pr, geom);

    //I'm temporarily hard coding in a constant B field to test
    Bsq = pow(200.,2.);

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

// void viscous_heating(){
//     // Need to find: u_ith  u_eth  u_nteth
//     double Tp = 

// }


#endif