#include "locmat.h"

int CellAvePorosity(const MeshInfo& mi, 
                    PhysProperty * pp,
                    basis& basis_,
                    const valarray<double>& gwf,
                    const vector<vertex>& gpf){

    double phi_f_hat = 0.0;
    double area = 0.0;

    for (unsigned int g=0; g<gwf.size(); g++){
        vertex mapped = GaussMapPointsFace(gpf[g],basis_.corners());
        double jac = abs(GaussJacobian(gpf[g],basis_.corners()));
        double gw = gwf[g];
        phi_f_hat += gw * jac * AssignPorosity(mapped,pp);
        area += gw * jac; 
    }

    phi_f_hat /= area;

    pp->phi_f_hat = phi_f_hat;

    return 0;
}

int AssignLocMat(const MeshInfo& mi,
                 BRMixed& br_,
                 basis& basis_,
                 LocMat * loc,
                 PhysProperty * pp,
                 const valarray<double>& gwe,
                 const valarray<double>& gpe,
                 const valarray<double>& gwf,
                 const vector<vertex>& gpf){

    // Cell average fluid porosity
    double phi_f_hat = pp->phi_f_hat;
    double phi_f = 0.0;
    double phi_s = 0.0;

    // Clear previous calculation
    loc->A.resize(12*12,0.0);
    loc->B.resize(12,0.0);
    loc->f.resize(12,0.0);

    std::fill(loc->A.begin(), loc->A.end(), 0.0);
    std::fill(loc->B.begin(), loc->B.end(), 0.0);
    std::fill(loc->f.begin(), loc->f.end(), 0.0);
    loc->C = 0.0;

    for (unsigned int g=0; g<gwf.size(); g++){
        // Calculate mapped gauss points and jacobian
        vertex mapped = GaussMapPointsFace(gpf[g],basis_.corners());
        double jac = abs(GaussJacobian(gpf[g],basis_.corners()));
        double gw = gwf[g];

        // Calculate point wise porosity ===================================
        phi_f = AssignPorosity(mapped, pp);  // Fluid porosity
        phi_s = AssignPorosity(phi_f);                 // Solid porosity

        std::array<std::array<double,4>, 12> brwork = 
                           br_.ComputeGradBRmixed(basis_, mapped);

        std::array<vertex, 12> brval = br_.ComputeBRmixed(basis_, mapped);

        vertex stokesforce = stokesForce(mapped,pp); 

        for (unsigned int j=0; j<12; j++){
                double div1 = brwork[j][0] + brwork[j][3];
            for (unsigned int i=0; i<12; i++){
                double A1 = brwork[j][0];
                double B1 = 0.5*(brwork[j][1] + brwork[j][2]);
                double C1 = brwork[j][3];

                double A2 = brwork[i][0];
                double B2 = 0.5*(brwork[i][1] + brwork[i][2]);
                double C2 = brwork[i][3];
 
                double div2 = brwork[i][0] + brwork[i][3];

                // Symmetrical formulation of A matrix
                loc->A[i+j*12] += 2*phi_s*gw*jac*
                                 (A1*A2+B1*B2*2+C1*C2 - (1.0/3.0)*div1*div2);
            }
            // With dimension version
            loc->B[j] += gw*jac*div1 * br_.Pressure();

            // Non dimensionalized version
            loc->f[j] += gw*jac* (stokesforce[0]*brval[j][0] + 
                                  stokesforce[1]*brval[j][1]);
        }

        // Non dimensionalized version
        loc->C += gw*jac*phi_f_hat/phi_s*
                  br_.Pressure()*br_.Pressure();

    }

    return 0;
}

int AssignLocMat(const MeshInfo& mi,
                 Hdivmixed& hdiv_,
                 basis& basis_,
                 LocMat * loc,
                 PhysProperty * pp,
                 const valarray<double>& gwe,
                 const valarray<double>& gpe,
                 const valarray<double>& gwf,
                 const vector<vertex>& gpf){

    // Cell average fluid porosity
    double phi_f_hat = pp->phi_f_hat;
    double phi_f = 0.0;
    double phi_s = 0.0;
    double theta = pp->theta;

    loc->A.resize(8*8,0.0);
    loc->B.resize(8,0.0);
    loc->f.resize(8,0.0);

    std::fill(loc->A.begin(), loc->A.end(), 0.0);
    std::fill(loc->B.begin(), loc->B.end(), 0.0);
    std::fill(loc->f.begin(), loc->f.end(), 0.0);
    loc->C = 0.0;

    for (unsigned int g=0; g<gwf.size(); g++){
        // Calculate mapped gauss points and jacobian
        vertex mapped = GaussMapPointsFace(gpf[g],basis_.corners());
        double jac = abs(GaussJacobian(gpf[g],basis_.corners()));
        double gw = gwf[g];

        // Calculate point wise porosity ===================================
        phi_f = AssignPorosity(mapped, pp);  // Fluid porosity
        phi_s = AssignPorosity(phi_f);

        std::array<vertex, 8> hdivwork = hdiv_.ComputeHdivmixed(basis_,mapped);

        vertex darcyforce = darcyForce(mapped,pp);

        for (unsigned int j=0; j<8; j++){
            for (unsigned int i=0; i<8; i++){

                // Non dimensionalized version
                loc->A[i+j*8] += gw*jac* 
                                 (hdivwork[j][0]*hdivwork[i][0] + 
                                  hdivwork[j][1]*hdivwork[i][1]);
            }
            // darctforce is set to be zero here
            loc->f[j] += gw*jac*(darcyforce[0]*hdivwork[j][0] + 
                                 darcyforce[1]*hdivwork[j][1]);
        }

        // Compaction matrix
        loc->C += gw*jac*1.0/phi_s*
                  hdiv_.Pressure()*hdiv_.Pressure();

    }

    // Define B matrix for the Darcy part
    // Compute with divergence theorem
    phi_f_hat = (phi_f_hat == 0.0 ? 1.0 : phi_f_hat);

    vertexSet corners = basis_.corners();

    for (int e =0; e<4; e++){

        vertexSet corner = {corners.at((e+3)%4),
                            corners.at(e)};
        double len = length(corner);
        for (int g=0; g<gpe.size(); g++){
            vertex mapped = GaussMapPointsEdge({gpe[g]},corner);
            std::array<vertex, 8>  hdivwork = hdiv_.ComputeHdivmixed(basis_,mapped);
            // Zeroth order constant pressure basis is always 1
            vertex nu = basis_.unitnormal(e);
            double phi_f_e = AssignPorosity(mapped, pp);
            for (int j=0; j<8; j++){
                // With dimension version
                loc->B[j] += len/2.0*gwe[g]*
                             pow(phi_f_hat,-0.5) * pow(phi_f_e, 1+theta) *
                             (hdivwork[j][0] * nu[0]+
                              hdivwork[j][1] * nu[1]);
            } 
        }
    }

    return 0;
}

int AssignLocMat(const MeshInfo& mi,
                 BRMixed& br_,
                 Hdivmixed& hdiv_,
                 basis& basis_,
                 PhysProperty * pp,
                 double * k,
                 const valarray<double>& gwf,
                 const vector<vertex>& gpf){

    *k = 0.0;

    double phi_f_hat = pp->phi_f_hat;
    double phi_f = 0.0;
    double phi_s = 0.0;

    for (unsigned int g=0; g<gwf.size(); g++){
        // Calculate mapped gauss points and jacobian
        vertex mapped = GaussMapPointsFace(gpf[g],basis_.corners());
        double jac = abs(GaussJacobian(gpf[g],basis_.corners()));
        double gw = gwf[g];

        // Calculate point wise porosity ===================================
        phi_f = AssignPorosity(mapped, pp);  // Fluid porosity
        phi_s = AssignPorosity(phi_f);

        *k -= gw*jac*pow(phi_f_hat,0.5)/phi_s * br_.Pressure() * 
                                                hdiv_.Pressure();
    }

    return 0;
}

