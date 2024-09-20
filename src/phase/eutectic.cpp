#include "eutectic.h"

EUTECTIC::phase::phase(){

    // Set up Clapeyron constant
    gamma_ = 10e-7;

    // Set up melting points under standard atmospheric pressure
	 // with dimension
	 Te0_ = 1560; //(K) 
	 T10_ = 2053; //(K) 

	 // dimensionless
	 // regardless of value of pressure
    TDe_ = 0;
	 TD1_ = 1;
	
	 // dimensionless latent heat 
    L_ = 0.3;

}

double EUTECTIC::phase::GetTD(const double& T, 
                              const double& P){

    // Get melting point with respect to current pressure P

    double Te = Te0_+gamma_*P;
    double T1 = T10_+gamma_*P;

    return (T-Te)/(T1-Te);
}

void EUTECTIC::phase::evalPhase(const double& HD,
                                const double& CD,
										  const double& P){

    switch(phaseSplit(HD,CD,P)){
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
            phi.mlt = HD/L_;
            phi.opx = CD - phi.mlt;
            phi.olv = 1-phi.mlt-phi.opx;
            TD      = gamma_*P;

            dTD_dCD = 0;
            dTD_dHD = 0;

        break;

        // Super eutectic two phase region
        case 4:
            TD      = ((HD+1) - sqrt(pow(HD+1,2)- 4*(HD-CD*L_)))/2;
            phi.opx = 0;
            phi.mlt = CD/(1-TD);
            phi.olv = 1-phi.opx-phi.mlt;

            dTD_dCD = -1./sqrt(pow(HD+1,2)- 4*(HD-CD*L_));
            dTD_dHD = 0.5 * (1 -1./sqrt(pow(HD+1,2)- 4*(HD-CD*L_)) * ((HD+1)-2));

        break;

        // Single phase super eutectic all melting region
        case 5:
            phi.opx = 0;
            phi.olv = 0;
            phi.mlt = 1;
            TD      = HD - L_;

            dTD_dCD = 0;
            dTD_dHD = 1;

        break;

        default:

            std::cout << "Invalid (H,C) pair." << " (" << HD << ", " << CD << ") " << std::endl;

        break;
    }
}

int EUTECTIC::phase::phaseSplit_(const double& HD, 
                                 const double& CD,
                                 const double& P){

    int region = 0;

    // Check valid pair of (C,H)
    assert(CD<1+std::numeric_limits<double>::epsilon());// "Opx composition beyond eutectic.\n");
    assert(CD>0-std::numeric_limits<double>::epsilon());// "Opx composition below zero.\n");

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

    return region;
}

int EUTECTIC::phase::phaseSplitTemp_(const double& TD,
                                     const double& CD,
                                     const double& phi2){

    int region = 0;


    return region;
}
