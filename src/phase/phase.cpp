#include "phase.h"

using namespace EUTECTIC;  

phase::phase(double X, double TD){
    // Set up phase physical properties eutectic using default physical values

    TDl_(X);
    Xhyb_(TD);

}

// Compute mass fraction
void FSolid1_(double X, double TD){
    fSolid1_ = (1-fSolid2_)*(TD<TDe) + 
               (1-fHybrid_)*((TD<TDL_ || TD==TDL_) && (TD>TDe || TD==TDe));
}

void FSolid2_(double X, double TD){
    fSolid2_ = X*(TD<TDe); 
}

void FHybrid_(double X, double TD){
    fHybrid_ = X/xHyb_*((TD<TDL_ || TD==TDL_) && (TD>TDe || TD==TDe)) + (TD>TDL_ || TD==TDL_);
}

// Compute volume fraction
void PhiSolid1_(double TD){
    phi::solid1 = (1-phi::solid2)*(TD<TDe) + 
                  (1-phi::hybrid)*((TD<TDL_ || TD==TDL_) && (TD>TDe || TD==TDe));
}

void PhiSolid2_(){
    phi::solid2 = fSolid2_/(hd::rho::ratio2 + (1-hd::rho::ratio2)*fSolid2_);
}

void PhiHybrid_(){
    phi::hybrid = fHybrid_/(hd::rho::ratio1 + (1-hd::rho::ratio1)*fHybrid_);
}

// dimensionless CD as function of TD and X
void CD_(){
    cd_ = phi::solid2*hd::rho::ratio2+
          phi::hybrid*hd::rho::ratio1*xHyb_;
}

void HD_(){
    hd_ = phi::solid1*hd::solid1 + 
          phi::solid2*hd::rho::ratio2*hd::solid2 + 
          phi::hybrid*hd::rho::ratio1*hd::hybrid; 
}

