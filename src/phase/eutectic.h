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

  T olv;
  T opx; 
  T mlt; 
};

class phase{

		  public:
            phase();
            ~phase() {};

            // Evaluate phase at certain pressure
            // Assign values to volumetric fractions
            // Evaluate non dimensionless values
            void evalPhase(const double& HD,
                           const double& CD);

            // Convert nondimensionlized variables to original variables
            void NonDimToDim(double pressure);

            // Volume fraction
            PhaseComp<double> phi;

            // Dimensionless temperature
            double TD;

        private:

            int phaseSplit(const double& H, 
                           const double& C);

            // Clapeyron constant
            // Relating perssure and melting point
            double gamma_;

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
};

}

#endif
