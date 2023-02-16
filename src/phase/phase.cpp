#include "phase.h"

using namespace EUTECTIC;  

void phase::SetPhase(double X, double TD){
    // Set up phase physical properties eutectic using default physical values
    TDl_(X);
    Xbri_(TD);

    MassFraction_(X,TD);

    VolumeFraction(TD);

    CD_();
    HD_(TD);
}

// Compute mass fraction
void phase::FIce_(double X, double TD){
    fIce_ = (1-fSal_)*(TD<TDe) + 
            (1-fBri_)*((TD<TDL_ || TD==TDL_) && (TD>TDe || TD==TDe));
}

void phase::FSal_(double X, double TD){
    fSal_ = X*(TD<TDe); 
}

void phase::FBri_(double X, double TD){
    fBri_ = X/xBri_*((TD<TDL_ || TD==TDL_) && (TD>TDe || TD==TDe)) + (TD>TDL_ || TD==TDL_);
}

// Compute volume fraction
void phase::PhiIce_(double TD){
    phi::ice = (1-phi::sal)*(TD<TDe) + 
               (1-phi::bri)*((TD<TDL_ || TD==TDL_) && (TD>TDe || TD==TDe));
}

void phase::PhiSal_(){
    phi::sal = fSal_/(hd::rho::si + (1-hd::rho::si)*fSal_);
}

void phase::PhiBri_(){
    phi::bri = fBri_/(hd::rho::bi + (1-hd::rho::bi)*fBri_);
}

// dimensionless CD as function of TD and X
void phase::CD_(){
    cd_ = phi::sal*hd::rho::si+
          phi::bri*hd::rho::bi*xBri_;
}

void phase::HD_(const double& TD){
    hd_ = phi::ice * hd::Ice(TD) + 
          phi::sal * hd::rho::si * hd::Sal(TD) + 
          phi::bri * hd::rho::bi * hd::Bri(TD,Ste); 
}

double invHX::HDe(const double& X, const double& Ste){
    return rho::bi/Ste*X/(rho::bi*Xe_ + (1-rho::bi)*X);
}

double invHX::HDl(const double& X, const double& Ste){
    return rho::bi*(1/Ste + cp::bi*(1-X/Xe_));
}

// Evaluation of current phase according to CD and HD values
evalPhase::evalPhase(const double& HD, const double& CD){
    HD_ = HD;
    CD_ = CD;

    if (HD < 1 && CD ==0){
        currentPhase_ = 1;
    } else if (HD<0 && CD>0){
        currentPhase_ = 2;
    } else if (HD>0 && HD<phase::invHC::HDe(CD,phase::Ste) && 
               CD>0 && CD < phase::invHC::CDb(HD,phase::Ste)){
        currentPhase_ = 3;
    } else if ((HD>phase::invHC::HDe(CD,phase::Ste) && HD<phase::invHC::HDl(CD,phase::Ste) && CD>0) || (HD>1 && HD<phase::invHC::HDl(CD,phase::Ste) && CD == 0)){
        currentPhase_ = 4;
    } else if (HD>phase::invHC::HDl(CD,phase::Ste)&&CD<phase::invHC::CD3l()) {
        currentPhase_ = 5;
    }

    SetPhase();

}

void evalPhase::SetPhase(){
    switch(currentPhase_){
        case 1: // No salt and no brine
            phase::phi::ice = 1;
            phase::phi::sal = 0;
            phase::phi::bri = 0;

            TD_ = HD_;

        case 2: // Two phase sub-solidus solid1 + solid2
            phase::phi::sal = CD_/phase::hd::rho::si;
            phase::phi::ice = 1 - phase::phi::sal;

            TD_ = HD_/(phase::phi::ice + phase::phi::sal) * 
                  phase::hd::rho::si * phase::hd::cp::si;

        case 3:
            phase::phi::bri = phase::Ste * HD_/phase::hd::rho::bi;
            phase::phi::sal = (CD_ - phase::Ste*phase::Xe*HD_)/phase::hd::rho::si;
            phase::phi::ice = 1 - phase::phi::bri - phase::phi::sal;

            TD_ = phase::TDe;

        case 4:
            if (CD_ == 0){
                TD_  = phase::invHC::field4T2(CD_,phase::Ste,HD_);
                phase::phi::bri = CD_/phase::hd::rho::bi*phase::Xe*(1-TD_);
            } else {
                TD_ = 1; // single phase limit of fieldion 4
                phase::phi::bri = (HD_-1)/(phase::hd::rho::bi/phase::Ste + 
                                           phase::hd::rho::bi*phase::hd::cp::bi-1);
            }

            // Volume fractions no salt
            phase::phi::ice = 1 - phase::phi::bri; 

        case 5:
            TD_ = phase::TDe + (HD_ - phase::hd::rho::bi/phase::Ste)/(phase::hd::rho::bi*
                                                                      phase::hd::cp::bi);
            phase::phi::bri = 1;

        default :
            std::cout << "Invalid phase..." << std::endl;
    }
}
