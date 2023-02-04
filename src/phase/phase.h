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

    // Dimensionless Specific enthalpy of the phases
    class hd : public rho, public cp{
        public:
            hd() {};
            ~hd() {};

            // Dimensionless specific enthalpy of the phases
            double solid1(const double& TD) const {return TD;};
            double solid2(const double& TD) const {return cp::ratio2*TD;};
            double hybrid(const double& TD, const double& Ste) const {return solid1(0.0) + 1/Ste + cp::ratio1*TD};

            // Bulk Enthalpy of system
            double HD(double TD, const phi& Phi) const {Phi.solid1*solid1(TD) +
                                                        Phi.solid2*rho::ratio2*solid2(TD) + 
                                                        Phi.hybrid1()*rho::ratio1*hybrid(TD)}; 

        private:
    }


    class phase : public cp{

        public: 
            phase() {};
            phase(double X, double TD);
            ~phase() {};

            double etutecticTemp  = 245; // Eutectic Temperature
            double multTemp1      = 273; // Multing Temperature of solid1 
            double multTemp2      = 400; // Multing Temperature of solid2
            double L              = 3.34e5; // Latent heat of water [J/kg]
 
            double DT = multTemp1 - etutecticTemp;
            double Ste = cp::solid1*DT/L; 

            double Xe;

        private:
            // Define liquidus and phase composition
            double TDe = 0;
            double TD1 = 1;

            void TDl_(double X) {TDL_ =  1-X/Xe;}; 

            void Xhyb_(double TD) {return Xe*(1-TD);}; 



            // Mass fractions
            void Fsolid1_(double X, double TD) {};
            void Fsolid2_(double X, double TD) {};
            void Fhybrid_(double X, double TD) {};

            double fsolid1_;
            double fsolid2_;
            double fhybrid_;

    }

   

}

#endif
