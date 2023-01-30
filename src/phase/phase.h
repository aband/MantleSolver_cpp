#ifndef PHASE_H_
#define PHASE_H_

namespace EUTECTIC{

    // Initial volume friction
    class phi{
        public:
            phi() {};
            ~phi() {}; 

            double solid1 = 0.3;
            double solid2 = 0.3;
            double hybrid1 = 0.2;
            double hybrid2 = 0.2;

    }

    // Density
    class rho{
        public:
            rho() {};
            ~rho() {}; 

            const double solid1 = 917; // kg/m^3
            const double solid2 = 1466; 
            const double hybrid = 1e3;

            // dimensionless density ratios
            const double ratio1 = hybrid/solid1;
            const double ratio2 = solid2/solid1;

    }

    // Heat capacity
    class cp{
        public:
            cp() {};
            ~cp() {};

            const double solid1 = 2000;
            const double solid2 = 920;
            const double hybrid = 4200;

            // dimensionless density ratios
            const double ratio1 = hybrid/solid1;
            const double ratio2 = solid2/solid1; 

    }

    // Thermal conductivity
    class kappa{
        public:
           kappa() {};
           ~kappa() {};

           void GetT(double T) {T_ = T;};

           const double Solid1() const {return 0.4685 + 488.12/T_;}; 
           const double Solid2() const {return 0.6  + T_*0;};
           const double Hybrid() const {return 0.56 + T_*0;};

           // dimensionless density ratios
           const double ratio1() const {return Hybrid()/Solid1();};
           const double ratio2() const {return Solid2()/Solid1();};
           const double ratio3() const {return Solid1()/Solid1();}; 

           const double sysD(const phi& Phi) const {Phi.solid1*ratio3() + 
                                                    Phi.solid2*ratio2() + 
                                                    Phi.hybrid1()*ratio1();};

        private:
           double T_;

    }

    // Dimensionless specific enthalpy of the phases
    class hd{
        public:
           hd() {};
           ~hd() {}; 

           
    }

    class phase : public phi, public rho, public kappa{
        public:
            phase() {};
            ~phase() {};

        private:
            double etutecticTemp_  = 245; // Eutectic Temperature
            double multTemp1_      = 273; // Multing Temperature of solid1 
            double multTemp2_      = 400; // Multing Temperature of solid2
            double L_              = 3.34e5; // Latent heat of water [J/kg]
 
            double DT_ = multTemp1_ - etutecticTemp_;
            double Ste = cp::solid1*DT_/L_; 

    }

}

#endif
