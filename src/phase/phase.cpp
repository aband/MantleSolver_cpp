#include "phase.h"

using namespace EUTECTIC;  

void phase::SetPhase(double X, double TD){
    // Set up phase physical properties eutectic using default physical values
    TDl_(X);
    Xbri_(TD);
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
