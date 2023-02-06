#ifndef PHASE_H_
#define PHASE_H_

namespace EUTECTIC{

    // Eutectic phase behavior exhibits 5 possible fieldions
    // 1) Pure           solid1
    // 2) sub-solidus    solid1 + solid2
    // 3) eutectic       solid1 + solid2 + brine
    // 4) super eutectic solid1 +           brine
    // 5) super liquidus                    brine

    // Initial volume friction
    class phi{
        public:
            phi() {};
            ~phi() {}; 

            double ice = 0.3;
            double sal = 0.3;
            double bri = 0.4;
            //double hybrid2 = 0.2;

    }

    // Density
    class rho{
        public:
            rho() {};
            ~rho() {}; 

            const double ice = 917; // kg/m^3
            const double sal = 1466; 
            const double bri = 1e3;

            // dimensionless density ratios
            const double bi = bri/ice;
            const double si = sal/ice;

    }

    // Heat capacity
    class cp{
        public:
            cp() {};
            ~cp() {};

            const double ice = 2000;
            const double sal = 920;
            const double bri = 4200;

            // dimensionless density ratios
            const double bi = bri/ice;
            const double si = sal/ice; 

    }

    // Thermal conductivity
    class kappa{
        public:
           kappa() {};
           ~kappa() {};

           void GetT(double T) {T_ = T;};

           const double Ice() const {return 0.4685 + 488.12/T_;}; 
           const double Sal() const {return 0.6  + T_*0;};
           const double Bri() const {return 0.56 + T_*0;};

           // dimensionless density ratios
           const double bi() const {return Bri()/Ice();};
           const double si() const {return Sal()/Ice();};
           const double ii() const {return Ice()/Ice();}; 

           const double sysD(const phi& Phi) const {Phi.ice*ii() + 
                                                    Phi.sal*si() + 
                                                    Phi.bri()*bi();};

        private:
           double T_;

    }

    // Dimensionless Specific enthalpy of the phases
    class hd : virtual public rho, public cp{
        public:
            hd() {};
            ~hd() {};

            // Dimensionless specific enthalpy of the phases
            double Ice(const double& TD) const {return TD;};
            double Sal(const double& TD) const {return cp::bi*TD;};
            double Bri(const double& TD, const double& Ste) const {return Ice(0.0) + 1/Ste + cp::bi*TD};

            // Bulk Enthalpy of system
            double HD(const double& TD, const phi& Phi) const {Phi.ice*Ice(TD) +
                                                               Phi.sal*rho::si * Sal(TD) + 
                                                               Phi.bri*rho::bi * Bri(TD)}; 

    }

    // Boundaries of the regions in HX-phase diagram
    class invHX : virtual public rho, virtual public cp{
        public : 
            invHX() {};
            ~invHX() {};

            double HDe(const double& X, const double& Ste, const double& Xe);

            double HDl(const double& X, const double& Ste, const double& Xe);

            double X1  = 0;
            double HD1 = 0;

            double X2s  = 0;
            double HD2s = 1;

            double X2l  = 0;
            double HD2l(const double& Ste) {return rho::bi*(1/Ste + cp::bi);};

            double X3s(const double& Xe) {return Xe;}
            double HD3s = 0;

            double X3l(const double& Xe) {return Xe;};
            double HD3l(const double& Ste) {return rho::bi/Ste;};

    }

    // Boundaries of the regions in HC-phase diagram
    class invHC : virtual public rho{
        public :
            invHC() {};
            ~invHC() {};

    }

    // Translating from matlab code to C++ code
	 // ice --- solid1
	 // salt --- solid2
	 // brine --- hybrid

    class phase : public phi, public hd{

        public: 
            phase() {};
            phase(double X, double TD);
            ~phase() {};

            double etutecticTemp  = 245; // Eutectic Temperature
            double multTemp1      = 273; // Multing Temperature of ice 
            double multTemp2      = 400; // Multing Temperature of sal
            double L              = 3.34e5; // Latent heat of water [J/kg]
 
            double DT = multTemp1 - etutecticTemp;
            double Ste = cp::ice*DT/L; 

            double Xe;

        private:
            // Define liquidus and phase composition
            double TDe = 0;
            double TD1 = 1;

            void TDl_(double X) {TDL_ =  1-X/Xe;}; 

            void Xbri_(double TD) {xBri_ = Xe*(1-TD);}; 

            double TDL_;
            double xBri_;

            // Mass fractions

            void MassFraction(double X, double TD) { FSal_(X,TD); FBri_(X,TD); FIce_(X,TD);};
            void FIce_(double X, double TD);
            void FSal_(double X, double TD);
            void FBri_(double X, double TD);

            double fIce_;
            double fSal_;
            double fBri_;

            // Volume fractions
            void VolumeFraction(double TD) { PhiSal_(); PhiBri_(); PhiIce_(TD);};
            void PhiIce_(double TD);
            void PhiSal_();
            void PhiBri_();

            // Dimensionless CD as function of TD and X
            void CD_();
 
            double cd_;

            // Dimensionless HD as function of TD and X
            void HD_();

            double hd_();
    }

   

}

#endif
