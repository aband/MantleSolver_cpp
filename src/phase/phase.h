#ifndef PHASE_H_
#define PHASE_H_

#include <math.h>
#include <iostream>

namespace EUTECTIC{

    // Eutectic phase behavior exhibits 5 possible fieldions
    // 1) Pure           solid1
    // 2) sub-solidus    solid1 + solid2
    // 3) eutectic       solid1 + solid2 + brine
    // 4) super eutectic solid1 +          brine
    // 5) super liquidus                   brine

    // Initial volume friction
    class phi{
        public:
            phi() {};
            ~phi() {}; 

            double ice = 0.3;
            double sal = 0.3;
            double bri = 0.4;
            //double hybrid2 = 0.2;

            double dicedxi = 0.0;
            double dsaldxi = 0.0;
            double dbridxi = 0.0;
    };

    // Density
    class rho{
        public:
            rho() {};
            ~rho() {}; 

            const double ice = 3e3; // kg/m^3
            const double sal = 3e3; 
            const double bri = 3e3;

            // dimensionless density ratios
            const double bi = bri/ice;
            const double si = sal/ice;

    };

    // Heat capacity
    class cp{
        public:
            cp() {};
            ~cp() {};

            const double ice = 1200;
            const double sal = 1200;
            const double bri = 1200;

            // dimensionless density ratios
            const double bi = bri/ice;
            const double si = sal/ice; 

    };

    // Thermal conductivity
    class kappa{
        public:
           kappa() {};
           ~kappa() {};

           void GetT(double T) {T_ = T;};

           //const double ice() const {return 0.4685 + 488.12/T_;}; 
           //const double sal() const {return 0.6  + T_*0;};
           //const double bri() const {return 0.56 + T_*0;};

           const double ice = 5.2;
           const double sal = 4.7;
           const double bri = 4.7;

           // dimensionless density ratios
           const double bi = bri/ice;
           const double si = sal/ice;
           const double ii = ice/ice;

           double sysD(const phi& Phi) const {return Phi.ice*ii + 
                                                     Phi.sal*si + 
                                                     Phi.bri*bi;};

        private:
           double T_;

    };

    // Dimensionless Specific enthalpy of the phases
    class hd : virtual public rho, virtual public cp{
        public:
            hd() {};
            ~hd() {};

            // Dimensionless specific enthalpy of the phases
            double Ice(const double& TD) const {return TD;};
            double Sal(const double& TD) const {return cp::bi*TD;};
            double Bri(const double& TD, const double& Ste) const {return Ice(0.0) + 1/Ste + cp::bi*TD;};

            // Bulk Enthalpy of system
            double HD(const double& TD, const phi& Phi, const double& Ste) const {return Phi.ice*Ice(TD) +
                                                                                         Phi.sal*rho::si * Sal(TD) + 
                                                                                         Phi.bri*rho::bi * Bri(TD, Ste);}; 

    };

    // Boundaries of the regions in HX-phase diagram
    class invHX : virtual public rho, virtual public cp{
        public : 
            invHX() {};
            ~invHX() {};

            void GetXe(const double& Xe) {Xe_ = Xe;};

            double HDe(const double& X, const double& Ste);

            double HDl(const double& X, const double& Ste);

            double X1  = 0;
            double HD1 = 0;

            double X2s  = 0;
            double HD2s = 1;

            double X2l  = 0;
            double HD2l(const double& Ste) {return rho::bi*(1/Ste + cp::bi);};

            double X3s() {return Xe_;}
            double HD3s = 0;

            double X3l() {return Xe_;};
            double HD3l(const double& Ste) {return rho::bi/Ste;};

            double field4T1(const double& X, const double& Ste, const double& HD) 
            {return -beta_(X,HD) + pow(discHX_(X,Ste,HD),0.5)/(2*alpha_(X));};

            double field4T2(const double& X, const double& Ste, const double& HD) 
            {return -beta_(X,HD) - pow(discHX_(X,Ste,HD),0.5)/(2*alpha_(X));};

            double disc(const double& X, const double& Ste, const double& HD) {return discHX_(X,Ste,HD);};

            double alpha(const double& X) {return alpha_(X);};
            double beta(const double& X, const double& HD) {return beta_(X, HD);};
            double gamma(const double& X, const double& Ste, const double& HD) 
            {return gamma_(X,Ste,HD);};

            private:
                double Xe_;
                // Solution to quadtatic in supra-eutectic region
                double a_(const double& X) {return rho::bi*Xe_/X + (1-rho::bi);};
                double b_(const double& X) {return -1*rho::bi*Xe_/X;};

                double alpha_(const double& X) {return b_(X);};
                double beta_(const double& X, const double& HD) 
                {return a_(X) + rho::bi*cp::bi - 1 - b_(X)*HD;};

                double gamma_(const double& X, const double& Ste, const double& HD)
                {return rho::bi/Ste-a_(X)*HD;};

                // Roots
                double discHX_(const double& X, const double& Ste, const double& HD)
                {return pow(beta_(X,HD),2) - 4*alpha_(X)*gamma_(X,Ste,HD);};
    };

    // Boundaries of the regions in HC-phase diagram
    class invHC : virtual public rho, virtual public cp{
        public :
            invHC() {};
            ~invHC() {};

            void GetXe(const double& Xe) {Xe_ = Xe;};

            double nu() {return (1-rho::si)*Xe_ + rho::si;};

            double HDe(const double& CD, const double& Ste) {return CD/(Xe_*Ste);};

            double HDl(const double& CD, const double& Ste) 
            {return rho::bi*(1/Ste + cp::bi * (1-CD/(rho::bi*Xe_)));};

            double HDb(const double& CD, const double& Ste)
            {return rho::bi*rho::si/((rho::bi*nu() - rho::si)*Ste *  (nu()*CD/(rho::si*Xe_)-1));};

            double CDb(const double& HD, const double& Ste)
            {return rho::si*Xe_/nu()*(1-Ste*HD/rho::bi) + Ste*Xe_*HD;};

            // Corners of phase fields in HC-diagram
            // 1
            double CD1 = 0;
            double HD1 = 0;
            // 2s
            double CD2s = 0;
            double HD2s = 1;
            // 2l
            double CD2l = 0;
            double HD2l(double invHXHD2l) {return invHXHD2l;};
            // 3s
            double CD3s() {return rho::si*Xe_/((1-rho::si)*Xe_+rho::si);};
            double HD3s = 0;
            // 3l
            double CD3l() {return rho::bi*Xe_;};
            double HD3l(double invHXHD3l) {return invHXHD3l;};

            double field4T1(const double& CD, const double& Ste, const double& HD)
            {return (-1*beta_(CD,HD) + pow(discHC_(CD,Ste,HD),0.5))/(2*alpha_());};

            double field4T2(const double& CD, const double& Ste, const double& HD)
            {return (-1*beta_(CD,HD) - pow(discHC_(CD,Ste,HD),0.5))/(2*alpha_());};

            double dfield4T2(const double& CD, const double& Ste, const double& HD)
            {return (-1*dbeta_(CD,HD) - 0.5*pow(discHC_(CD,Ste,HD),-0.5)*ddiscHC_(CD,Ste,HD))/(2*alpha_());};

            double disc(const double& CD, const double& Ste, const double& HD)
            {return discHC_(CD,Ste,HD);};

            double alpha(const double& CD, const double& HD)
            {return alpha_();};

            double beta(const double& CD, const double& HD)
            {return beta_(CD,HD);};

            double gamm(const double& CD, const double& Ste, const double& HD)
            {return gamma_(CD, Ste, HD);}

            private:
                double Xe_;
                double alpha_() {return rho::bi*Xe_;};
                double beta_(const double& CD, const double& HD)
                {return (1-rho::bi*cp::bi)*CD - (1+HD)*rho::bi*Xe_;};
                double dbeta_(const double& CD, const double& HD)
                {return (1-rho::bi*cp::bi);};

                double gamma_(const double& CD, const double& Ste, const double& HD)
                {return rho::bi*Xe_*HD - rho::bi/Ste*CD;};
                double dgamma_(const double& CD, const double& Ste, const double& HD)
                {return  -1*rho::bi/Ste;};

                double discHC_(const double& CD, const double& Ste, const double& HD)
                {return pow(beta_(CD,HD),2) - 4*alpha_() * gamma_(CD,Ste,HD);};
                double ddiscHC_(const double& CD, const double& Ste, const double& HD)
                {return 2*beta_(CD,HD)*dbeta_(CD,HD) - 4*alpha_()*dgamma_(CD, Ste, HD);};
    };

    // Translating from matlab code to C++ code
	 // ice --- solid1
	 // salt --- solid2
	 // brine --- hybrid

    class phase : public phi, public hd, public invHX, public invHC{

        public: 
            phase() {};
            ~phase() {};

            void SetPhase(double X, double TD);

            double eutecticTemp = 1227; // Eutectic Temperature
            double multTemp1    = 1350; // Multing Temperature of ice 
            double multTemp2    = 400; // Multing Temperature of sal
            double L            = 4e5; // Latent heat of water [J/kg]
 
            double DT = multTemp1 - eutecticTemp;

            //void UpdateSte() {Ste = cp::ice*DT/L;};
            double Ste = cp::ice*DT/L; 

            double Xe = 0.7;

            // Define liquidus and phase composition
            double TDe = 0;
            double TD1 = 1;

        private:
            void TDl_(double X) {TDL_ =  1-X/Xe;}; 

            void Xbri_(double TD) {xBri_ = Xe*(1-TD);}; 

            double TDL_;
            double xBri_;

            // Mass fractions

            void MassFraction_(double X, double TD) { FSal_(X,TD); FBri_(X,TD); FIce_(X,TD);};
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
            void HD_(const double& TD);

            double hd_;
    };

    class evalPhase : public phase , public kappa{
        public:
        // Identify different fieldions of HX phase diagram
        evalPhase(const double& HD, const double& CD);
        ~evalPhase();

        void EvalPhase(const double& HD, const double& CD);
        // Evaluate volume fraction according to different regimes
        void SetPhase();

        double getD1s(){return hd::rho::ice;};
        double getD2s(){return hd::rho::sal;};
        double getD1f(){return hd::rho::ice;};
        double getD2f(){return hd::rho::sal;};

        // Get three volume fraction
        double getPhi1(){return phase::phi::ice;};
        double getPhi2(){return phase::phi::sal;};
        double getPhif(){return phase::phi::bri;};

        double getDPhi1(){return phase::phi::dicedxi;};
        double getDPhi2(){return phase::phi::dsaldxi;};
        double getDPhif(){return phase::phi::dbridxi;};

        void ViewPhysics();

        // Functions used for testing
        void ViewPhase();

        private:
        int currentPhase_ = 0;
        double CD_;
        double HD_;
        double TD_;
    };

}

#endif
