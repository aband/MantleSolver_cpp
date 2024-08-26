#ifndef EUTECTIC_H_
#define EUTECTIC_H_

// Final version of eutectic 
#include <math.h>
#include <iostream>
#include <cmath>
#include <limits>

namespace EUTECTIC{

template <typename T>
struct PhaseComp{

  // A template struct holding information regarding three components 
  // in the eutectic phase package
  // For example, volumetric values for three components

  T solid1;
  T solid2; 
  T fluid; 
};

class phase{

		  public:
            phase();
            ~phase() {};

            // Evaluate phase at certain pressure
				// Assign values to volumetric fractions
				// Evaluate non dimensionless values
            void evalPhase();

            // Convert nondimensionlized variables to original variables
            void NonDimToDim(double pressure);

            // Volume fraction
				PhaseComp<double> phi;

        private:

            int phaseSplit(double H, double C);

            // Clapeyron constant
				// Relating perssure and melting point
            double gamma_;

            // Melting temperatures under standard atmospheric pressure
				// with dimension
			   double te0_;
			   double t10_;	

            // Melting temperatures under standard atmospheric pressure
				// dimensionless
            double Te0_;
				double T10_;

            // Melting temperatures (dimensionless)
            double Te_;
            double T1_;

            // Latent heat
				// dimensionless
            double L_;

};

}

#endif
