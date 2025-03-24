#ifndef EUTECTIC_RESCALED_H_
#define EUTECTIC_RESCALED_H_

// New eutectic phase package
// Modify with proper scale again

#include <math.h>
#include <iostream>
#include <cmath>
#include <limits>
#include <cassert>

namespace EUTECTIC{

    struct PhaseComp{

        // A template struct holding information regarding three components 
        // in the eutectic phase package
        // For example, volumetric values for three components
        // (Attention, T does not represent temperature).

        // volumetric fraction
        double phi1;
        double phi2; 
        double phil;

        // mass fraction
        double cl;

        // Phase region
        int region;

        // Melting temperature
        double TDp;

        // Derivative
        double dTD_dCD;
        double dTD_dHD; 
    };

    class phase{

        public:
            phase();
            ~phase() {};

            double GetTDp(const double& T, 
                          const double& P) const;      // Compute pressure corrected temperature points

            double GetStaticP(const double& zD,
                              const double& l0) const; // Compute static pressure with dimension 

            int evalPhase(const double& HD,
                          const double& CD,
                          const double& P);

            int phaseSplit(const double& inHD,
                           const double& inCD,
                           const double& P);

            PhaseComp pc;

            int printInfo() const;

            double Tm0;   // Standard melting point
            double Te0;   // Standard eutectic point
            double nu;    // Clapeyron constant
            double gamma; // inverse Clapeyron constant
            double dT;    // Temperature difference between melting and eutectic temperature

            double L;     // Latent heat
            double LD;    // Dimensionless Latent heat
            double g;     // Gravity accelaration

            double cp;    // Specific heat capacity
            double TDm0;  // Dimensionless standard melting temperature
            double TDe0;  // Dimensionless standard eutectic temperature
            double rho;   // Density
            double rhor;  // Density difference

            double Xe;    // Eutectic liquid composition

            double mus;
            double mul;
       
            double k0;
            double invk0;

            double l0;
            double u0;
            double p0;
            double t0;
            double alpha0;
    };

}

#endif
