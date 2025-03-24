#include "eutectic_rescaled.h"

EUTECTIC::phase::phase(){

    nu    = 6.5*pow(10,6);
    gamma = 1.0/nu;

    Tm0   = 2000;
    Te0   = 1480;
    dT    = Tm0 - Te0;

    L     = 5*pow(10,5);
    cp    = 1200;

    LD    = L/cp/dT;
    TDm0  = Tm0/dT;
    TDe0  = Te0/dT;

    rho   = 3000;
    rhor  = 500;
    Xe    = 0.25;

    mus   = 1e19;
    mul   = 1;

    k0    = 1e-8;
    invk0 = 1e8;

    g     = 10;

    l0    = pow(mus*k0/mul,0.5);
    u0    = k0*rhor*g/mul;
    p0    = rhor*g*l0;
    t0    = l0/u0;
    
    alpha0 = 3e-5;
}

double EUTECTIC::phase::GetTDp(const double& TD, 
                               const double& P) const{
    return TD + gamma * P/dT;
}

double EUTECTIC::phase::GetStaticP(const double& zD,
                                   const double& l0) const{
    return rho*g*zD*l0;
}

// Taking dimensionless enthalpy, dimensionless composition and
// Pressure with dimension as input
// Be careful!
int EUTECTIC::phase::evalPhase(const double& HD,
                               const double& CD,
                               const double& P){

    int region = phaseSplit(HD, CD, P);

    double Tep = GetTDp(TDe0, P);
    double Tmp = GetTDp(TDm0, P);

    switch(region){
        
        case 1:
            pc.phi1 = 1;
            pc.phi2 = 0;
            pc.phil = 0;

            if (HD > Tmp) {
                pc.TDp = Tmp;
                pc.dTD_dHD = 0.0;
            }else {
                pc.TDp = HD;
                pc.dTD_dHD = 1.0;
            } 

            pc.dTD_dCD = 0.0;

            pc.cl = 0.0;
        break;

        case 2:
            pc.phi1 = 1-CD;
            pc.phi2 = CD;
            pc.phil = 0;

            pc.TDp  = HD; 
            pc.dTD_dCD = 0.0;
            pc.dTD_dHD = 1.0;

            pc.cl = 0.0;
        break;

        case 3:
            pc.phil = (HD - Tep)/ LD;
            pc.phi2 = CD - pc.phil * Xe;
            pc.phi1 = 1 - pc.phil - pc.phi2;

            pc.TDp  = Tep;

            pc.dTD_dCD = 0;
            pc.dTD_dHD = 0;

            pc.cl = Xe;

        break;

        case 4:
            pc.TDp  = 0.5*(HD + Tmp - sqrt(pow(HD - Tmp,2) + 4*CD*LD/Xe));
            pc.phi2 = 0.0;
            pc.phil = CD/(Xe*(Tmp - pc.TDp)); 
            pc.phi1 = 1 - pc.phi2 - pc.phil;
          
            pc.dTD_dCD = -LD/Xe * 1.0/sqrt( pow(HD - Tmp,2) + 4*CD*LD/Xe );
            pc.dTD_dHD = 0.5* (1 - (HD - Tmp)/sqrt( pow(HD - Tmp,2) + 4*CD*LD/Xe ) );

            pc.cl = Xe*(Tmp - pc.TDp);
        break;

        case 5:
            pc.phi1 = 0.0;
            pc.phi2 = 0.0; 
            pc.phil = 1.0;

            pc.TDp  = HD - LD;

            pc.dTD_dCD = 0.0;
            pc.dTD_dHD = 1.0; 
            pc.cl = CD;
        break;

        default:

            std::cout << "Invalid (H,C) pair." << " (" << HD << ", " << CD << ") " << std::endl;

        break;
    }

    return 1;
}

int EUTECTIC::phase::phaseSplit(const double& inHD,
                                const double& inCD,
                                const double& P){

    int region = 0;

    double HD = inHD;
    double CD = inCD;
    // Alter dimensionless HD and CD 
    if (inCD > Xe){
        CD = Xe;
    }

    if (inCD < 0){
        CD = 0;
    }

    // Compute current eutectic and melting points
    double Tep = GetTDp(TDe0, P);
    double Tmp = GetTDp(TDm0, P);

    // Two lines separating phase regions
    double lineb = Tep + LD * CD/Xe;                         // line separating sub and super eutectic regions
    double linec = Tmp - CD/Xe + LD; // line separating super eutectic and all meltiing region

    if (CD < std::numeric_limits<double>::epsilon() && 
        HD < Tmp + LD){
        // Pure component 1
        region = 1;
    }

    if (HD <= Tep && CD > 0){
        // Sub-eutectic solid phase
        region = 2;
    } else if (HD < lineb + std::numeric_limits<double>::epsilon() &&
               HD > Tep && CD > 0){
        // Eutectic region
        region = 3;
    } else if (HD > lineb && 
               HD < linec + std::numeric_limits<double>::epsilon() && CD > 0){
        // Super-eutectic region
        region = 4;
    } else if (HD > linec && CD > 0){
        // All melting 
        region = 5;
    }

    pc.region = region;

    return region;
}

int EUTECTIC::phase::printInfo() const{

    std::cout << "The eutectic phase package has been defined ..." << std::endl;
    std::cout << "Dimensionless melting temperature  : " << TDm0 << std::endl;
    std::cout << "Dimensionless eutectic temperature : " << TDe0 << std::endl;
    std::cout << "Dimensionless latent heat          : " << LD     << std::endl;
    std::cout << "Characteristic length              : " << l0     << std::endl;
    std::cout << "Characteristic velocity            : " << u0     << std::endl;
    std::cout << "Characteristic pressure            : " << p0     << std::endl;
    std::cout << "Characteristic time                : " << t0     << std::endl;

    return 1;
}
