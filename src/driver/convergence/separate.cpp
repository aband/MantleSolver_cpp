int Driver::CellAvePorosity_convergence(const indice& gcell){

    const valarray<double>& gwf = GaussWeightsFace;
    const vector<vertex>&   gpf = GaussPointsFace;

    double phi_f_hat = 0.0;
    double area = 0.0;

    //std::string cellLoc = location(mi, gcell);

    for (unsigned int g=0; g<gwf.size(); g++){
        vertex mapped = GaussMapPointsFace(gpf[g],basis_->corners());

        // Get HD and CD from reconstruction at this gaussian point
        double phif = 0.0;

        // Test =================================================
        phif = AssignPorosity(mapped, myPhase->pp); 
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

int Driver::AssignLocMatStokes_convergence(const indice& gcell,
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

    //phi_f_hat = (phi_f_hat == 0.0 ? 1.0 : phi_f_hat);

    for (unsigned int g=0; g<gwf.size(); g++){
        // Calculate mapped gauss points and jacobian
        vertex mapped = GaussMapPointsFace(gpf[g],basis_->corners());
        double jac = abs(GaussJacobian(gpf[g],basis_->corners()));
        double gw = gwf[g];

        // Reconstruction of point wise value of HD and CD
        // Test ==================================================================
        phi_f = AssignPorosity(mapped, myPhase->pp);
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
//        loc->C += gw*jac*phi_f/phi_s*
//                  br_->Pressure()*br_->Pressure();
    }

    return 1;
}


