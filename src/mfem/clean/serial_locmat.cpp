#include "serial_solver.h"

// Porosity will not be treated or averaged here, 
// but in separate porosity function instead

static int clearLocMat(int size,
                       LocMat& loc){

    loc.A.resize(size*size, 0.0);
    loc.B.resize(size, 0.0);
    loc.f.resize(size, 0.0);
    
    std::fill(loc.A.begin(), loc.A.end(), 0.0);
    std::fill(loc.B.begin(), loc.B.end(), 0.0);
    std::fill(loc.f.begin(), loc.f.end(), 0.0);
    loc.C = 0.0;

    return 1;
}

int DarcyStokes::AssignLocMatStokes(const MeshInfo& mi,
                                    LocMat& loc,
                                    double theta, 
                                    const poroSet& poro,
                                    PhysProperty * pp){

    clearLocMat(12, loc);
   
    // copy gaussian quadrature points
    const valarray<double>& gwe = GaussWeightsEdge;
    const valarray<double>& gpe = GaussPointsEdge;
    const valarray<double>& gwf = GaussWeightsFace;
    const vector<vertex>&   gpf = GaussPointsFace;

    for (int g=0; g<gwf.size(); g++){
        // Calculate mapped gauss points and jacobian
        vertex mapped = GaussMapPointsFace(gpf[g],basis_.corners());
        double jac = abs(GaussJacobian(gpf[g],basis_.corners()));
        double gw = gwf[g];

        std::array<std::array<double,4>, 12> brwork = 
                           br_.ComputeGradBRmixed(basis_, mapped);

        std::array<vertex, 12> brval = br_.ComputeBRmixed(basis_, mapped);

        vertex stokesforce = stokesForce(mapped,pp); 

        double phif = poro.cellporo.at(g); 
        double phis = 1-phif;

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
                loc.A.at(i+j*12) += 2*phis*gw*jac*
                                      (A1*A2+B1*B2*2+C1*C2 - (1.0/3.0)*div1*div2);
            }
            // With dimension version
            loc.B.at(j) += gw*jac*div1 * br_.Pressure();

            // Non dimensionalized version
            // Attention, porosity has been multiplied to right hand side force term
            loc.f.at(j) += gw*jac* phis*(stokesforce[0]*brval[j][0] + 
                                         stokesforce[1]*brval[j][1]);

        }

        // Non dimensionalized version
//        loc->C += gw*jac*phi_f_hat/phis*
//                  br_->Pressure()*br_->Pressure();
        loc.C += gw*jac*phif/phis*
                 br_.Pressure()*br_.Pressure();
    }

    return 1;
}

int DarcyStokes::AssignLocMatDarcy(const MeshInfo& mi,
                                   LocMat& loc,
                                   double theta,
                                   const poroSet& poro,
                                   PhysProperty * pp){

    // copy gaussian quadrature points
    const valarray<double>& gwe = GaussWeightsEdge;
    const valarray<double>& gpe = GaussPointsEdge;
    const valarray<double>& gwf = GaussWeightsFace;
    const vector<vertex>&   gpf = GaussPointsFace;

    clearLocMat(8, loc);  

    for (unsigned int g=0; g<gwf.size(); g++){
        // Calculate mapped gauss points and jacobian
        vertex mapped = GaussMapPointsFace(gpf[g],basis_.corners());
        double jac = abs(GaussJacobian(gpf[g],basis_.corners()));
        double gw = gwf[g];

        double phif = poro.cellporo.at(g); 
        double phis = 1-phif;

        std::array<vertex, 8> hdivwork = hdiv_.ComputeHdivmixed(basis_,mapped);

        vertex darcyforce = darcyForce(mapped,pp);

        for (unsigned int j=0; j<8; j++){
            for (unsigned int i=0; i<8; i++){

                // Non dimensionalized version
                loc.A.at(i+j*8) += gw*jac* 
                                  (hdivwork[j][0]*hdivwork[i][0] + 
                                   hdivwork[j][1]*hdivwork[i][1]);
            }
            // darctforce is set to be zero here
            loc.f.at(j) += gw*jac*(darcyforce[0]*hdivwork[j][0] + 
                                   darcyforce[1]*hdivwork[j][1]);
        }

        // Compaction matrix
        double scaletmp = 1.0;
        if (poro.aveporo > 1e-15){
            scaletmp = phif/poro.aveporo;
        }
        loc.C += gw*jac*scaletmp/phis*
                 hdiv_.Pressure()*hdiv_.Pressure();
    }

    double phif_hat = (poro.aveporo == 0.0 ? 1.0 : poro.aveporo);

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

//cout <<"Edge : " << e << "  porosity : " << g << " is : " <<  poro.edgeporo.at(e*gpe.size()+g) << endl;
            for (int j=0; j<8; j++){
                // With dimension version
                loc.B[j] += len/2.0*gwe[g]*
                            pow(phif_hat,-0.5) * pow(poro.edgeporo.at(e*gpe.size()+g), 1+theta) *
                            (hdivwork[j][0] * nu[0]+
                             hdivwork[j][1] * nu[1]);
            } 
        }
    }

    return 1;
}

int DarcyStokes::AssignLocMatCouple(const MeshInfo& mi,
                                    double& k,
                                    double theta,
                                    const poroSet& poro,
                                    PhysProperty * pp){

    const valarray<double>& gwf = GaussWeightsFace;
    const vector<vertex>&   gpf = GaussPointsFace;

    k = 0.0;

    double phi_f = 0.0;
    double phi_s = 0.0;

    for (unsigned int g=0; g<gwf.size(); g++){
        // Calculate mapped gauss points and jacobian
        vertex mapped = GaussMapPointsFace(gpf[g],basis_.corners());
        double jac = abs(GaussJacobian(gpf[g],basis_.corners()));
        double gw = gwf[g];

        // Reconstruction of point wise value of HD and CD

        phi_f = poro.cellporo.at(g);

        phi_s = 1.0-phi_f;

        //k -= gw*jac*pow(phi_f_hat,0.5)/phi_s * br_->Pressure() * 
        //                                       hdiv_->Pressure();

        double scaletmp = 0.0;
        if (poro.aveporo > 1e-16){
            scaletmp = phi_f/sqrt(poro.aveporo);
        }

        k -= gw*jac*scaletmp/phi_s * br_.Pressure() * 
                                     hdiv_.Pressure();
    }

    return 1;
}
