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
        double c1l;
        double c2l;

        // Phase region
        int region;

        // Melting temperature
        double Tm_p;

        // Derivative
        T dTD_dCD;
        T dTD_dHD; 
    };

    class phase{

        public:
            phase();
            ~phase() {};

            double GetTDp(const double& T, 
                          const double& P) const;      // Compute pressure corrected temperature points

            double GetStaticP(const double& zD,
                              const double& l0) const; // Compute static pressure 
 
        private:
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
    };

}

#endif
