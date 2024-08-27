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

            // Get corresponding dimensionless temperature with respect to each pressure
            double GetTD(const double& T, 
                         const double& P);

            // Split phase regions according to values of dimensionless enthalpy and 
            // dimensionless composition.
            int phaseSplit(const double& HD,
                           const double& CD) {return phaseSplit_(HD, CD);};

            // Split phase regions according to value of dimensionless temperature and
            // dimensionless composition
            // Need phi2 for exact determination
            int phaseSplitTemp(const double& TD,
                               const double& CD,
                               const double& phi2) {return phaseSplitTemp_(TD,CD,phi2);};

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

            int phaseSplit_(const double& HD, 
                            const double& CD);

            int phaseSplitTemp_(const double& TD,
                                const double& CD,
                                const double& phi2);

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
};

}

#endif
