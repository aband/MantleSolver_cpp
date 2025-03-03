#include "eutectic.h"

EUTECTIC::phase::phase(){

    // Set up Clapeyron constant
    gamma_ = 10e-7;

    // Set up melting points under standard atmospheric pressure
	 // with dimension
	 Te0_ = 1480; //(K) 
	 T10_ = 2053; //(K) 

	 // dimensionless
	 // regardless of value of pressure
    TDe_ = 0;
	 TD1_ = 1;
	
	 // dimensionless latent heat 
    L_ = 0.3;

    // Uniform density
	 // No need to differentiate densities of two minerals
    rho_ = 3000; //kg/m^3

    g_   = 10;

}

double EUTECTIC::phase::GetTe(const double& P){

    // Get melting point with respect to current pressure P

    return TD*(T10_-Te0_) + Te0_ + gamma_*P;
}

double EUTECTIC::phase::GetScaledLithoP(const double& z){

    // Return lithostatic pressure scaled with (T1-Te)
    // input z of unit [m]
    //return rho_*g_*z/(T10_ - Te0_);

    return 5.23 * z;
}

double EUTECTIC::phase::GetDepth(const double& y, const double& l0){

    return y*(-1)*l0*0.6;
}

double EUTECTIC::phase::Getef(){

    return phi.mlt*(TD + L_);
}

double EUTECTIC::phase::Getcf(const double& CD){
    return CD - phi.opx;

}

void EUTECTIC::phase::evalPhase(const double& HD,
                                const double& CD,
                                const double& P){

    // Pressure corrected melting temperature
    double Tm = 1+ gamma_*P;
    phi.Tm_p = gamma_*P;

    switch(phaseSplit(HD,CD,P)){
        // Single phase solidus
        case 1:
            phi.olv = 1;
            phi.opx = 0;
            phi.mlt = 0;
            TD      = HD;

            phi.dTD_dCD = 0;
            phi.dTD_dHD = 1;

        break;
 
        // Two phase solidus
        case 2:
            phi.opx = CD;
            phi.olv = 1-CD;
            phi.mlt = 0;
            TD      = HD;

            phi.dTD_dCD = 0;
            phi.dTD_dHD = 1;

        break;

        // Three phase eutectic
        case 3:
            phi.mlt = (HD-gamma_*P)/L_;
            phi.opx = CD - phi.mlt;
            phi.olv = 1-phi.mlt-phi.opx;
            TD      = gamma_*P;

            phi.dTD_dCD = 0;
            phi.dTD_dHD = 0;

        break;

        // Super eutectic two phase region
        case 4:
            TD      = ((HD+Tm) - sqrt(pow(HD+Tm,2)- 4*(Tm*HD-CD*L_)))/2;
            phi.opx = 0;
            phi.mlt = CD/(Tm-TD);
            phi.olv = 1-phi.opx-phi.mlt;

            phi.dTD_dCD = -L_/sqrt(pow(HD+Tm,2)- 4*(Tm*HD-CD*L_));
            phi.dTD_dHD = 0.5 * (1 + 0.5/sqrt(pow(HD+Tm,2)- 4*(Tm*HD-CD*L_)) * (2*HD - 4*Tm) );

        break;

        // Single phase super eutectic all melting region
        case 5:
            phi.opx = 0;
            phi.olv = 0;
            phi.mlt = 1;
            TD      = HD - L_;

            phi.dTD_dCD = 0;
            phi.dTD_dHD = 1;

        break;

        default:

            std::cout << "Invalid (H,C) pair." << " (" << HD << ", " << CD << ") " << std::endl;

        break;
    }
}

void EUTECTIC::phase::evalPhase(const double& HD, 
                                const double& CD,
                                const double& P,
                                const double& rho_f,
                                const double& rho_s){

    double Tm = 1+ gamma_*P;

    double a = 0;
    double b = 0;
    double c = 0;

    switch(phaseSplit_(HD,CD,P,rho_f,rho_s)){
        // Single phase solidus
        case 1:
            phi.olv = 1;
            phi.opx = 0;
            phi.mlt = 0;
            TD      = HD;

            dTD_dCD = 0;
            dTD_dHD = 1;

        break;
 
        // Two phase solidus
        case 2:
            phi.opx = CD;
            phi.olv = 1-CD;
            phi.mlt = 0;
            TD      = HD;

            dTD_dCD = 0;
            dTD_dHD = 1;

        break;

        // Three phase eutectic
        case 3:
            phi.mlt = rho_s/rho_f*(HD-gamma_*P)/L_;
            phi.opx = CD - rho_f/rho_s*phi.mlt;
            phi.olv = 1-phi.mlt-phi.opx;
            TD      = gamma_*P;

            dTD_dCD = 0;
            dTD_dHD = 0;

        break;

        // Super eutectic two phase region
        case 4:
            a = 1;
            b = HD + Tm - (1-rho_f/rho_s)*rho_s/rho_f*CD;
            c = Tm*HD - CD*L_;

				TD      = (b - sqrt(b*b - 4*a*c))/2;
            phi.opx = 0;
            phi.mlt = rho_s/rho_f*CD/(Tm-TD);
            phi.olv = 1-phi.opx-phi.mlt;

            dTD_dCD = -1./sqrt(pow(HD+Tm,2)- 4*(Tm*HD-CD*L_));
// It is WRONG!!!! 
//            dTD_dHD = 0.5 * (1 -1./sqrt(pow(HD+1,2)- 4*(HD-CD*L_)) * ((HD+1)-2));

        break;

        // Single phase super eutectic all melting region
        case 5:
            phi.opx = 0;
            phi.olv = 0;
            phi.mlt = 1;
            TD      = rho_s/rho_f*HD - L_;

            dTD_dCD = 0;
            dTD_dHD = 1;

        break;

        default:

            std::cout << "Invalid (H,C) pair." << " (" << HD << ", " << CD << ") " << std::endl;

        break;
    }
   
}

int EUTECTIC::phase::phaseSplit_(const double& inHD, 
                                 const double& inCD,
                                 const double& P){

    int region = 0;

    // Check valid pair of (C,H)
//    assert(CD<1+std::numeric_limits<double>::epsilon());// "Opx composition beyond eutectic.\n");
//    assert(CD>0-std::numeric_limits<double>::epsilon());// "Opx composition below zero.\n");

    double HD = inHD;
    double CD = inCD;
    if(inCD > 1) {
        CD = 1;
    }

    if (inCD < 0){
        CD = 0;
    }

    // Two lines separating phase regions
    double lineb = L_* CD + gamma_ * P;
    double linec = 1+L_-CD + gamma_ * P;

    if (CD < std::numeric_limits<double>::epsilon() && 
        HD < 1+std::numeric_limits<double>::epsilon() + gamma_*P){
        // Single phase sub-solidus region
        region = 1;

    } else if (HD < std::numeric_limits<double>::epsilon() + gamma_*P){
        // Two phase sub-soidus region
        region = 2;

    } else if (HD < lineb + std::numeric_limits<double>::epsilon() && 
               HD > 0 + gamma_*P){
        // Three phase eutectic region
        region = 3;
    } else if (HD > lineb && 
               HD < linec + std::numeric_limits<double>::epsilon()){
        // Two phase eutectic region
        region = 4;
    } else if (HD > linec){
        // All melting
        region = 5;
    }

    phi.region = region;

    return region;
}

int EUTECTIC::phase::phaseSplit_(const double& HD, 
                                 const double& CD,
                                 const double& P,
                                 const double& rho_f,
                                 const double& rho_s){

    int region = 0;

    double coef = rho_f/rho_s;

    // Rescale CD with coef inverse

    region = phaseSplit_(HD, CD/coef, P);

    return region;
}
