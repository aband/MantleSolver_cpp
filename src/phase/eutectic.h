#ifndef EUTECTIC_H_
#define EUTECTIC_H_

// Final version of eutectic 
#include <math.h>
#include <iostream>
#include <cmath>
#include <limits>
#include <cassert>

namespace EUTECTIC{

template <typename T>
struct PhaseComp{

  // A template struct holding information regarding three components 
  // in the eutectic phase package
  // For example, volumetric values for three components
  // (Attention, T does not represent temperature).

  T olv;
  T opx; 
  T mlt;

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

            // Get dimentionalized temperature 
            double GetTe(const double& P);

            // Get lithostatic pressure 
				// Take depth as input and return lithostatic pressure
            double GetScaledLithoP(const double& z);

            double GetDepth(const double& y, 
                            const double& l0);

            double Getef();
            double Getcf(const double& CD);

            // Split phase regions according to values of dimensionless enthalpy and 
            // dimensionless composition.
            int phaseSplit(const double& HD,
                           const double& CD,
                           const double& P) {return phaseSplit_(HD, CD, P);};

            int phaseSplit(const double& HD,
                           const double& CD) {return phaseSplit_(HD, CD, 0);};

            // Evaluate phase at certain pressure
            // Assign values to volumetric fractions
            // Evaluate non dimensionless values
            void evalPhase(const double& HD,
                           const double& CD,
                           const double& P);

            void evalPhase(const double& HD,
                           const double& CD) {evalPhase(HD,CD,0);};

            void evalPhase(const double& HD, 
                           const double& CD,
                           const double& P,
                           const double& rho_f,
                           const double& rho_s);

            // Convert nondimensionlized variables to original variables
            void NonDimToDim(double pressure);

            // Volume fraction
            PhaseComp<double> phi;

            //PhaseComp<double> c;

            //PhaseComp<double> e;

            // Dimensionless temperature
            double TD;

            // Derivative of dimensionless temperature against
            // concentration C and enthalpy.
            double dTD_dCD;

            double dTD_dHD;

        private:

            int phaseSplit_(const double& HD, 
                            const double& CD,
                            const double& P);
 
            // Correction separating liquid and solid densities
            int phaseSplit_(const double& HD,
                            const double& CD,
                            const double& P,
                            const double& rho_f,
                            const double& rho_s);

            // Clapeyron constant
            // Relating perssure and melting point
            double gamma_; // K*pa^-1

            // Melting temperatures under standard atmospheric pressure
            // with dimension
            double Te0_;
            double T10_;	

            // Melting temperatures (dimensionless)
            double TDe_;
            double TD1_;

            // Latent heat
            // dimensionless
            double L_;

            // Density
            double rho_;

            // Gravitional acceleration
            double g_;
};

}

#endif
