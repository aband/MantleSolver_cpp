#ifndef SIMPLEPHASE_H_
#define SIMPLEPHASE_H_

// A simplified phase code
#include <math.h>
#include <iostream>

template <typename T>
struct PhaseComp{

  // A template struct holding information regarding three components 
  // in the eutectic phase package

  T solid1;
  T solid2; 
  T fluid; 
};

class phaseState{
    public:
        phaseState();
        ~phaseState();

        int GetVolumeFrac(PhaseComp<double> Phi) {phi.solid1 = Phi.solid1;
                                                  phi.solid2 = Phi.solid2;
                                                  phi.fluid  = Phi.fluid};

        int GetVolumeFrac(double solid1, double solid2, double fluid) 
                         {phi.solid1 = solid1;
                          phi.solid2 = solid2;
                          phi.fluid  = fluid};

        PhaseComp<double> phi; // Volume fraction

        int EvalPhase();

    private:

        // Constants
        double L_;  // Latent heat in melting 
        double Xe_; // Eutectic composition of solid 2
        double T1_; // Melting temperature for solid 1
        double Te_; // Eutectic temperature

        // Reference physical attributes used for nondimensionalization 
        double DT_;
        double hc_;
        double Hc_;
        double Cc_;

        PhaseComp<double> rho_;   // Density
        PhaseComp<double> cp_;    // Specific enthalpy
        PhaseComp<double> kappa_; // Thermal conductivity

        double ste_;  // Stefan number (ratio of sensible heat to latent heat)

        // Points and lines divides phase diagram
        double CD1_ ;
        double HD1_ ;
        double CD2s_;
        double HD2s_;
        double CD2l_;
        double HD2l_;
        double CD3s_;
        double HD3s_;
        double CD3l_;
        double HD3l_;
};

#endif
