#include "driver.h"

/**
 * Customized transport function
 * Transporting phi with customized melting rate directly
 */

// Two sided
int Driver::computeEffVel_case(const vector<vertex>& gaussp,
                               const vertexSet& edgep,
                               const indice& gcellin, const indice& gcellout,
                               const Tensor<weights>& allwgts, double ** lphi,
                               vector<vertex>& vel){
 
    // Extract velocity on these given gauss points
    vector<vertex> vel_relative = 
    ExtractVelocity(&sresult_->vel_darcy, &sresult_->g_darcy,
                    refArrayDarcyEssen_,mi,
                    gaussp, gcellin,*hdiv_,*basis_,{1});
    
    vector<vertex> vel_stokes = 
    ExtractVelocity(&sresult_->vel_stokes, &sresult_->g_stokes,
                    refArrayStokesEssen_,mi,
                    gaussp, gcellin,*br_,*basis_,{1});

    for (int g=0; g<gaussp.size(); g++){
        double phi_mean = 0.0;

        double phiin = advection.eval(gaussp.at(g), ml, location(mi,gcellin), 
                       allwgts({gcellin[0], gcellin[1]}), gcellin, lphi);

        double phiout = advection.eval(gaussp.at(g), ml, location(mi,gcellout), 
                        allwgts({gcellout[0], gcellout[1]}), gcellout, lphi);

        if (phiin == 0.0 && phiout == 0.0){
            phi_mean = 0.0;
        } else {
            phi_mean  = harmonic_mean(phiin.at(g) ,phiout.at(g));
        }

        vel.at(g) = phi_mean*(vel_relative.at(g) + vel_stokes.at(g));

    }

    return 1;
}

// One sided
int Driver::computeEffVel_case(const vector<vertex>& gaussp,
                               const vertexSet& edgep,
                               const indice& gcell,
                               const Tensor<weights>& allwgts, double ** lphi,
                               vector<vertex>& vel){
 
    // Extract velocity on these given gauss points
    vector<vertex> vel_relative = 
    ExtractVelocity(&sresult_->vel_darcy, &sresult_->g_darcy,
                    refArrayDarcyEssen_,mi,
                    gaussp, gcellin,*hdiv_,*basis_,{1});
    
    vector<vertex> vel_stokes = 
    ExtractVelocity(&sresult_->vel_stokes, &sresult_->g_stokes,
                    refArrayStokesEssen_,mi,
                    gaussp, gcellin,*br_,*basis_,{1});

    for (int g=0; g<gaussp.size(); g++){

        double phi = advection.eval(gaussp.at(g), ml, location(mi,gcell), 
                     allwgts({gcell[0], gcell[1]}), gcell, lphi);

        vel.at(g) = phi*(vel_relative.at(g) + vel_stokes.at(g));

    }

    return 1;
}

int Driver::updateEdgeFlux_case(Tensor<double>& vertedge, Tensor<double>& horiedge,
                                const Tensor<weights>& allwgts, double ** lphi){

    Tensor_zero(vertedge);
    Tensor_zero(horiedge);

    const valarray<double>& gwe = GaussWeightsEdge;
    const valarray<double>& gpe = GaussPointsEdge;

    std::vector<vertex> gaussp;
    gaussp.resize(gpe.size());

    // Effective velocity should be used for transport of concentration
    vector<vertex> effvel; effvel.resize(gaussp.size());

    vector<double> uin;  uin.resize(gaussp.size());
    vector<double> uout; uout.resize(gaussp.size());
 
    for (int j=0; j<mi.MPIglobalCellSize[1]; j++){
    for (int i=0; i<mi.MPIglobalCellSize[0]; i++){

        double flux = 0.0;

        indice gcell {i,j} = 0.0;
        indice cellout;

        effvel.clear(); effvel.resize(gaussp.size()); 

        // Extract corners with respect to given global indice
        vertexSet corners = extractCorners(mi, gcell); 

        // =============================================================
        // Get horizontal edge
        vertexSet hori {corners.at(0), corners.at(1)};

        // Extract velocity on this edge 
        for (int g=0; g<gpe.size(); g++){
            gaussp.at(g) = GaussMapPointsEdge({gpe[g]},hori);
        }   

        if (j==0){
           
           computeEffVel(gaussp, hori, gcell, allwgts, lphi, effvel);

           flux = edgefluxintegral(hori, CDbottom ,effvel);

        } else {

           computeEffVel(gaussp, hori, gcellin, cellout, allwgts, lphi, effvel);

           flux = edgefluxintegral(mi, gcell, cellout, hori, allwgts, effvel, ml, advection, lphi);

        }        

        horiedge({i,j}) = flux;

        // 1D problem
        vertedge({i,j}) = 0.0;   

    }}

    for (int i=0; i<mi.MPIglobalCellSize[0]; i++){

        // Regarded as outside cell
        indice gcell {i, mi.MPIglobalCellSize[1]-1};
        vertexSet corners = extractCorners(mi, gcell);
        vertexSet hori    = {corners.at(3), corners.at(2)};

        effvel.clear(); effvel.resize(gaussp.size()); 

        // Extract velocity on this edge 
        for (int g=0; g<gpe.size(); g++){
            gaussp.at(g) = GaussMapPointsEdge({gpe[g]},hori);
        }   

        double flux = edgefluxintegral(mi, gcell, hori, allwgts, effvel, ml, advection, lphi);

        horiedge({i, mi.MPIglobalCellSize[1]}) = flux; 
    }

    return 1;
}

// ===== Assembly ====
int Driver::CellAvePorosity_case(const indice& gcell, 
                                 const Tensor<weights>& allwgts,
                                 double ** lphi){

    const valarray<double>& gwf = GaussWeightsFace;
    const vector<vertex>&   gpf = GaussPointsFace;

    double phi_f_hat = 0.0;
    double area = 0.0;

    //std::string cellLoc = location(mi, gcell);

    for (unsigned int g=0; g<gwf.size(); g++){
        vertex mapped = GaussMapPointsFace(gpf[g],basis_->corners());

        // Get HD and CD from reconstruction at this gaussian point
        double phif = advection.eval(mapped, ml, location(mi,gcell), allwgts({gcell[0], gcell[1]}), gcell, lphi);

        // Test =================================================

        //phif = AssignPorosity(mapped, myPhase->pp); 

        // ======================================================

        double jac = abs(GaussJacobian(gpf[g],basis_->corners()));
        double gw = gwf[g];
        phi_f_hat += gw * jac * phif;
        area += gw * jac; 
    }

    phi_f_hat /= area;

    myPhase->pp->phi_f_hat = phi_f_hat;

    return 1;
}

int Driver::AssignLocMatStokes_case(const indice& gcell,
                                    const Tensor<weights>& allwgts,
                                    double ** lphi,
                                    LocMat * loc){

    // copy gaussian quadrature points
    const valarray<double>& gwe = GaussWeightsEdge;
    const valarray<double>& gpe = GaussPointsEdge;
    const valarray<double>& gwf = GaussWeightsFace;
    const vector<vertex>&   gpf = GaussPointsFace;

    // Cell average fluid porosity
    double phi_f_hat = myPhase->pp->phi_f_hat;
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
        vertex mapped = GaussMapPointsFace(gpf[g],basis_->corners());
        double jac = abs(GaussJacobian(gpf[g],basis_->corners()));
        double gw = gwf[g];

        // Reconstruction of point wise value of HD and CD
        double phi_f = advection.eval(mapped, ml, location(mi,gcell), allwgts({gcell[0], gcell[1]}), gcell, lphi);

        // Test ==================================================================

        //phi_f = AssignPorosity(mapped, myPhase->pp);

        // =======================================================================

        phi_s = AssignPorosity(phi_f);       // Solid porosity

        std::array<std::array<double,4>, 12> brwork = 
                           br_->ComputeGradBRmixed(*basis_, mapped);

        std::array<vertex, 12> brval = br_->ComputeBRmixed(*basis_, mapped);

        vertex stokesforce = stokesForce(mapped,myPhase->pp); 

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
            loc->B[j] += gw*jac*div1 * br_->Pressure();

            // Non dimensionalized version
            // Attention, porosity has been multiplied to right hand side force term
            loc->f[j] += gw*jac* phi_f*(stokesforce[0]*brval[j][0] + 
                                        stokesforce[1]*brval[j][1]);

        }

        // Non dimensionalized version
        loc->C += gw*jac*phi_f_hat/phi_s*
                  br_->Pressure()*br_->Pressure();

    }

    return 1;
}

int Driver::AssignLocMatDarcy_case(const indice& gcell,
                                   const Tensor<weights>& allwgts,
                                   double ** lphi,
                                   LocMat * loc){
                     
    // copy gaussian quadrature points
    const valarray<double>& gwe = GaussWeightsEdge;
    const valarray<double>& gpe = GaussPointsEdge;
    const valarray<double>& gwf = GaussWeightsFace;
    const vector<vertex>&   gpf = GaussPointsFace;

    // Cell average fluid porosity
    double phi_f_hat = myPhase->pp->phi_f_hat;
    double phi_f = 0.0;
    double phi_s = 0.0;
    double theta = myPhase->pp->theta;

    loc->A.resize(8*8,0.0);
    loc->B.resize(8,0.0);
    loc->f.resize(8,0.0);

    std::fill(loc->A.begin(), loc->A.end(), 0.0);
    std::fill(loc->B.begin(), loc->B.end(), 0.0);
    std::fill(loc->f.begin(), loc->f.end(), 0.0);
    loc->C = 0.0;

    for (unsigned int g=0; g<gwf.size(); g++){
        // Calculate mapped gauss points and jacobian
        vertex mapped = GaussMapPointsFace(gpf[g],basis_->corners());
        double jac = abs(GaussJacobian(gpf[g],basis_->corners()));
        double gw = gwf[g];

        // Reconstruction of point wise value of HD and CD
        double phi_f = advection.eval(mapped, ml, location(mi, gcell), allwgts({gcell[0], gcell[1]}), gcell, lphi);

        // Test ==================================================================

        //phi_f = AssignPorosity(mapped, myPhase->pp);

		  // =======================================================================

        phi_s = AssignPorosity(phi_f);

        std::array<vertex, 8> hdivwork = hdiv_->ComputeHdivmixed(*basis_,mapped);

        vertex darcyforce = darcyForce(mapped,myPhase->pp);

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
                  hdiv_->Pressure()*hdiv_->Pressure();

    }

    // Define B matrix for the Darcy part
    // Compute with divergence theorem
    phi_f_hat = (phi_f_hat == 0.0 ? 1.0 : phi_f_hat);

    vertexSet corners = basis_->corners();

    for (int e =0; e<4; e++){

        vertexSet corner = {corners.at((e+3)%4),
                            corners.at(e)};
        double len = length(corner);
        for (int g=0; g<gpe.size(); g++){
            vertex mapped = GaussMapPointsEdge({gpe[g]},corner);
            std::array<vertex, 8>  hdivwork = hdiv_->ComputeHdivmixed(*basis_,mapped);
            // Zeroth order constant pressure basis is always 1
            vertex nu = basis_->unitnormal(e);

            // Reconstruction of point wise value of HD and CD
            double phi_f_e = advection.eval(mapped, ml, location(mi, gcell), allwgts({gcell[0], gcell[1]}), gcell, lphi);

            // Testing =================================================

            //phi_f_e = AssignPorosity(mapped, myPhase->pp);

            // =========================================================

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

int Driver::AssignLocMatCouple_case(const indice& gcell,
                                    const Tensor<weights>& allwgts,
                                    double ** lphi,
                                    double& k){

    const valarray<double>& gwf = GaussWeightsFace;
    const vector<vertex>&   gpf = GaussPointsFace;

    k = 0.0;

    double phi_f_hat = myPhase->pp->phi_f_hat;
    double phi_f = 0.0;
    double phi_s = 0.0;

    for (unsigned int g=0; g<gwf.size(); g++){
        // Calculate mapped gauss points and jacobian
        vertex mapped = GaussMapPointsFace(gpf[g],basis_->corners());
        double jac = abs(GaussJacobian(gpf[g],basis_->corners()));
        double gw = gwf[g];

        // Reconstruction of point wise value of HD and CD
        double phi_f = advection.eval(mapped, ml, location(mi, gcell), allwgts({gcell[0], gcell[1]}), gcell, lphi);

        // Test ============================================================

        //phi_f = AssignPorosity(mapped, myPhase->pp);

        // =================================================================

        phi_s = AssignPorosity(phi_f);

        k -= gw*jac*pow(phi_f_hat,0.5)/phi_s * br_->Pressure() * 
                                               hdiv_->Pressure();
    }

    return 1;
}

int Driver::PrepareTransport_case(double (*func)(const valarray<double>& point,
                                                 const vector<double>& param) ){

    PetscCall(DMCreateGlobalVector(dmu, &globalCD));

    SimpleInitialValue(dmMesh, dmu, &globalmesh, &globalCD, {H_,0.0}, func);

    // Initialization of multi level weno and corresponding usage
    ml = multilevel(); 

    ml.addLevel("(1,3)", {1,3}, mi);
    ml.addLevel("(1,2)", {1,2}, mi);

    // Area scale
    h0 = sqrt((L_*H_)/
         (double)(mi.MPIglobalCellSize[0]*mi.MPIglobalCellSize[1]));

    advection = mluse();

    unordered_map<std::string, vector<indice>> method;
    method.insert(std::make_pair<std::string, vector<indice>>("(1,3)", { {0,-1} }));
    method.insert(std::make_pair<std::string, vector<indice>>("(1,2)", { {0,-1} , {0,0} }));

    advection.setmethod("all", method);
    advection.setbias("all");

    CDbottom = func({0.0,-1*H_}, {H_,0.0});

    return 1;
}

int Driver::RK_case(double dt, double Tmax, int maxIter, double tolUzawa){

    int mark = 1 + start;

    int Nt = (int)(Tmax/dt);

    for (int t=0; t<Nt; t++) {

        Vec localphi;
        PetscCall(DMGetLocalVector(dmu, &localphi));

        Vec fluxphi;
        PetscCall(VecDuplicate(globalCD, &fluxphi));
 
        double ** lphi;
        double ** lfphi;

        PetscCall(DMGlobalToLocalBegin(dmu, globalCD, INSERT_VALUES, localphi));
        PetscCall(DMGlobalToLocalEnd(dmu, globalCD, INSERT_VALUES, localphi));

        PetscCall(DMDAVecGetArray(dmu, localphi, &lphi););

        PetscCall(DMDAVecGetArray(dmu, fluxphi, &lfphi));




    }

    return 1;
}
