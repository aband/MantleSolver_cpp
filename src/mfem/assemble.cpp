#include "assemble.h"

void AssignLocMatrix(const MeshInfo& mi,
                     basis& basis_,
                     Hdivmixed& hdiv_,
                     BRMixed& br_,
                     LocMatrix * locmatrix,
                     PhysProperty * physproperty,
                     const valarray<double>& gwe,
                     const valarray<double>& gpe,
                     const valarray<double>& gwf,
                     const vector<vertex>& gpf){

    // Temporary physics parameters
    double theta  = physproperty->theta;
    double mu_s   = physproperty->mu_s;
    double mu_f   = physproperty->mu_f;
    double inv_k0 = 1;
    double rho_r  = physproperty->rho_f/physproperty->rho_s;
    double gx     = physproperty->gx;
    double gy     = physproperty->gy;

    double phi_f = 1.0;
    double phi_s = 0.0;

    // Cell average fluid porosity
    double phi_f_hat = 0.0;
    double area = 0.0;

    // Perpare with area and cell averaged porosity
    for (unsigned int g=0; g<gwf.size(); g++){
        vertex mapped = GaussMapPointsFace(gpf[g],basis_.corners());
        double jac = abs(GaussJacobian(gpf[g],basis_.corners()));
        double gw = gwf[g];
        phi_f_hat += gw * jac * AssignPorosity(mapped,(*physproperty).l);
        area += gw * jac; 
    }

    phi_f_hat /= area;

    // Assign values to local matrix
    // Zeros out all local values first
    for (unsigned int j=0; j<8; j++){(*locmatrix).bd[j] = 0.0;
                                     (*locmatrix).sourcedarcy[j] = 0.0;
        for (unsigned int i=0; i<8; i++){(*locmatrix).ad[j*8+i] = 0.0;}}

    for (unsigned int j=0; j<12; j++){(*locmatrix).bs[j] = 0.0;
                                      (*locmatrix).sourcestokes[j] = 0.0;
        for (unsigned int i=0; i<12; i++){(*locmatrix).as[j*12+i] = 0.0;}}

    (*locmatrix).cs = 0.0;
    (*locmatrix).cd = 0.0;
    (*locmatrix).k  = 0.0;

    for (unsigned int g=0; g<gwf.size(); g++){
        vertex mapped = GaussMapPointsFace(gpf[g],basis_.corners());
        double jac = abs(GaussJacobian(gpf[g],basis_.corners()));
        double gw = gwf[g];

        std::array<std::array<double,4>, 12> brwork = 
                           br_.ComputeGradBRmixed(basis_, mapped);

        std::array<vertex, 12> brval = br_.ComputeBRmixed(basis_, mapped);

        phi_f = AssignPorosity(mapped, (*physproperty).l);
        phi_s = 1 - phi_f;

        vertex stokesforce = stokesForce(mapped); 

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

                // Symmetrical formulation
                //(*locmatrix).as[i+j*12] += gw*jac* 2*mu_s*phi_s * (A1*A2+B1*B2*2+C1*C2);
                //(*locmatrix).as[i+j*12] += gw*jac*2*(A1*A2+B1*B2*2+C1*C2 - (1.0/3.0)*div1*div2);

                // Nonsymmetrical formulation
                (*locmatrix).as[i+j*12] += gw*jac*(brwork[j][0]*brwork[i][0] + 
                                                   brwork[j][1]*brwork[i][1] + 
                                                   brwork[j][2]*brwork[i][2] + 
                                                   brwork[j][3]*brwork[i][3]);
            }

            (*locmatrix).bs[j] += gw*jac*div1 * 1;

            //(*locmatrix).sourcestokes[j] += gw*jac*(1-phi_f)*rho_r*(stokesforce[0]*brwork[j][0] + 
            //                                                        stokesforce[1]*brwork[j][1]);

            (*locmatrix).sourcestokes[j] += gw*jac*(stokesforce[0]*brval[j][0] + 
                                                    stokesforce[1]*brval[j][1]);
        }

        std::array<vertex, 8> hdivwork = hdiv_.ComputeHdivmixed(basis_,mapped);

        vertex darcyforce = darcyForce(mapped);

        for (unsigned int j=0; j<8; j++){
            for (unsigned int i=0; i<8; i++){
                (*locmatrix).ad[i+j*8] += gw*jac*mu_f*inv_k0* 
                                (hdivwork[j][0]*hdivwork[i][0] + 
                                 hdivwork[j][1]*hdivwork[i][1]);
            }
            (*locmatrix).sourcedarcy[j] += gw*jac*(darcyforce[0]*hdivwork[j][0] + 
                                                   darcyforce[1]*hdivwork[j][1]);
        }

        (*locmatrix).cd += gw*jac*1.0/(mu_s*(1-phi_f))*1*1;

        (*locmatrix).cs += gw*jac*phi_f_hat/(mu_s*(1-phi_f))*1*1;

        (*locmatrix).k += gw*jac*pow(phi_f_hat,0.5) * 1*1;
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
            for (int j=0; j<8; j++){
                (*locmatrix).bd[j] += len/2.0*gwe[g]*
                                      pow(phi_f_hat,-0.5) * pow(phi_f, 1+theta) *
                                     (hdivwork[j][0] * nu[0]+
                                      hdivwork[j][1] * nu[1]);
            } 
        }
    }
}

PetscErrorCode SerialMatrixAssembleBlock(const MeshInfo& mi,
                                         basis& basis_,
                                         Hdivmixed& hdiv_,
                                         BRMixed& br_,
                                         PhysProperty * physproperty,
                                         System * system){
 
    PetscErrorCode    ierr;
    PetscFunctionBeginUser;

    const valarray<double>& gwe = GaussWeightsEdge;
    const valarray<double>& gpe = GaussPointsEdge;
    const valarray<double>& gwf = GaussWeightsFace;
    const vector<vertex>& gpf = GaussPointsFace;

    PetscCall(MatCreate(PETSC_COMM_WORLD,&(*system).Ad));
    PetscCall(MatCreate(PETSC_COMM_WORLD,&(*system).As));
    PetscCall(MatCreate(PETSC_COMM_WORLD,&(*system).Bs));
    PetscCall(MatCreate(PETSC_COMM_WORLD,&(*system).Bd));
    PetscCall(MatCreate(PETSC_COMM_WORLD,&(*system).Cs));
    PetscCall(MatCreate(PETSC_COMM_WORLD,&(*system).Cd));
    PetscCall(MatCreate(PETSC_COMM_WORLD,&(*system).K));

    PetscCall(VecCreate(PETSC_COMM_WORLD,&system->sourceStokes));
    PetscCall(VecCreate(PETSC_COMM_WORLD,&system->sourceDarcy));

    PetscCall(VecSetSizes(system->sourceStokes, PETSC_DECIDE, br_.getDOF()));
    PetscCall(VecSetSizes(system->sourceDarcy, PETSC_DECIDE, hdiv_.getDOF()));

    PetscCall(VecSetUp(system->sourceStokes));
    PetscCall(VecSetUp(system->sourceDarcy));

    double *sourcestokes;
    double *sourcedarcy;

    PetscCall(VecGetArray(system->sourceStokes, &sourcestokes));
    PetscCall(VecGetArray(system->sourceDarcy, &sourcedarcy));

    // Set zeros to right hand side vector
    for (unsigned int k=0; k<br_.getDOF(); k++){sourcestokes[k] = 0.0;}

    for (unsigned int k=0; k<hdiv_.getDOF(); k++){sourcedarcy[k] = 0.0;}

    int totalElem = mi.MPIglobalCellSize[0] * mi.MPIglobalCellSize[1];

    PetscCall(MatSetSizes((*system).Ad,PETSC_DECIDE,PETSC_DECIDE,
                                       hdiv_.getDOF(),hdiv_.getDOF()));
    PetscCall(MatSetSizes((*system).As,PETSC_DECIDE,PETSC_DECIDE,
                                       br_.getDOF(),br_.getDOF()));

    PetscCall(MatSetSizes((*system).Bd,PETSC_DECIDE,PETSC_DECIDE,
                                       hdiv_.getDOF(),totalElem));
    PetscCall(MatSetSizes((*system).Bs,PETSC_DECIDE,PETSC_DECIDE,
                                       br_.getDOF(),totalElem));

    PetscCall(MatSetSizes((*system).Cd,PETSC_DECIDE,PETSC_DECIDE,
                                       totalElem,totalElem));
    PetscCall(MatSetSizes((*system).Cs,PETSC_DECIDE,PETSC_DECIDE,
                                       totalElem,totalElem));

    PetscCall(MatSetSizes((*system).K,PETSC_DECIDE,PETSC_DECIDE,
                                       totalElem,totalElem));

    PetscCall(MatSetType((*system).Ad,MATMPIAIJ));
    PetscCall(MatSetType((*system).As,MATMPIAIJ));
    PetscCall(MatSetType((*system).Bd,MATMPIAIJ));
    PetscCall(MatSetType((*system).Bs,MATMPIAIJ));
    PetscCall(MatSetType((*system).Cd,MATMPIAIJ));
    PetscCall(MatSetType((*system).Cs,MATMPIAIJ));
    PetscCall(MatSetType((*system).K,MATMPIAIJ));
 
    ierr = MatSetUp((*system).As);CHKERRQ(ierr);
    ierr = MatSetUp((*system).Bs);CHKERRQ(ierr);
    ierr = MatSetUp((*system).Bd);CHKERRQ(ierr);
    ierr = MatSetUp((*system).Cs);CHKERRQ(ierr);
    ierr = MatSetUp((*system).Cd);CHKERRQ(ierr);
    ierr = MatSetUp((*system).K);CHKERRQ(ierr);
    ierr = MatSetUp((*system).Ad);CHKERRQ(ierr);

    //LocMatrix * locmatrix = (LocMatrix *)malloc(sizeof(LocMatrix));

    LocMatrix * locmatrix = new LocMatrix;

    int checkSizeM, checkSizeN;
    MatGetSize((*system).Ad, &checkSizeM, &checkSizeN);

    for (unsigned int n=0; n<totalElem; n++){
        indice globalElem = Bend(mi,n);

        basis_.GetCorners(mi,globalElem); 

        std::array<int, 8> tmp = hdiv_.LocalToGlobal(mi,globalElem);

        // Compute local matrix
        AssignLocMatrix(mi,basis_,hdiv_,br_,locmatrix,physproperty,gwe,gpe,gwf,gpf);

        const int ISDarcy[8] = 
           {tmp[0],tmp[1],tmp[2],tmp[3],tmp[4],tmp[5],tmp[6],tmp[7]};

        for (unsigned int k=0; k<8; k++){
            sourcedarcy[ISDarcy[k]] += (*locmatrix).sourcedarcy[k];
        }

        const double Bdv[8] = {(*locmatrix).bd[0],(*locmatrix).bd[1],
                               (*locmatrix).bd[2],(*locmatrix).bd[3],
                               (*locmatrix).bd[4],(*locmatrix).bd[5],
                               (*locmatrix).bd[6],(*locmatrix).bd[7]}; 

        const int idxm = n;

        PetscCall(MatSetValuesBlocked((*system).Bd, 8, ISDarcy, 1, &idxm, Bdv, ADD_VALUES));

        for (unsigned int l=0; l<8; l++){
            const double Adv[8] = {(*locmatrix).ad[0+l*8],(*locmatrix).ad[1+l*8],
                                   (*locmatrix).ad[2+l*8],(*locmatrix).ad[3+l*8],
                                   (*locmatrix).ad[4+l*8],(*locmatrix).ad[5+l*8],
                                   (*locmatrix).ad[6+l*8],(*locmatrix).ad[7+l*8]};

            const int idxm = ISDarcy[l];
            PetscCall(MatSetValuesBlocked((*system).Ad,1,&idxm,8,ISDarcy,Adv,ADD_VALUES));
        }

        PetscCall(MatSetValue((*system).Cd,n,n,(*locmatrix).cd,ADD_VALUES));

        // Assemble Stokes part

        const std::array<int, 12> tmp2 = br_.LocalToGlobal(mi,globalElem);

        const int ISStokes[12] = {tmp2[0],tmp2[1],tmp2[2],tmp2[3],
                                  tmp2[4],tmp2[5],tmp2[6],tmp2[7],
                                  tmp2[8],tmp2[9],tmp2[10],tmp2[11]};

        for (unsigned int k=0; k<12; k++){
            sourcestokes[ISStokes[k]] += (*locmatrix).sourcestokes[k];
        }

        const double Bsv[12] = {(*locmatrix).bs[0],(*locmatrix).bs[1],
                                (*locmatrix).bs[2],(*locmatrix).bs[3],
                                (*locmatrix).bs[4],(*locmatrix).bs[5],
                                (*locmatrix).bs[6],(*locmatrix).bs[7],
                                (*locmatrix).bs[8],(*locmatrix).bs[9],
                                (*locmatrix).bs[10],(*locmatrix).bs[11]};

        const int idxms = n;
        PetscCall(MatSetValuesBlocked((*system).Bs,12,ISStokes,1,&idxms,Bsv,ADD_VALUES));

        for (unsigned int l=0; l<12; l++){
            const double Asv[12] = {(*locmatrix).as[0+l*12], (*locmatrix).as[1+l*12],
                                    (*locmatrix).as[2+l*12], (*locmatrix).as[3+l*12],
                                    (*locmatrix).as[4+l*12], (*locmatrix).as[5+l*12],
                                    (*locmatrix).as[6+l*12], (*locmatrix).as[7+l*12],
                                    (*locmatrix).as[8+l*12], (*locmatrix).as[9+l*12],
                                    (*locmatrix).as[10+l*12], (*locmatrix).as[11+l*12]};

            const int idxm = ISStokes[l];
            PetscCall(MatSetValuesBlocked((*system).As, 1, &idxm, 12, ISStokes, Asv, ADD_VALUES));
        }

        PetscCall(MatSetValue((*system).Cs, n, n, (*locmatrix).cs, ADD_VALUES));

        PetscCall(MatSetValue((*system).K, n, n, (*locmatrix).k, ADD_VALUES));

    }

    PetscCall(VecRestoreArray(system->sourceStokes,&sourcestokes));
    PetscCall(VecRestoreArray(system->sourceDarcy,&sourcedarcy));

    // Assemble all block matrix
    PetscCall(MatAssemblyBegin((*system).As,MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd((*system).As,MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyBegin((*system).Ad,MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd((*system).Ad,MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyBegin((*system).Bs,MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd((*system).Bs,MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyBegin((*system).Bd,MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd((*system).Bd,MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyBegin((*system).Cs,MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd((*system).Cs,MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyBegin((*system).Cd,MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd((*system).Cd,MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyBegin((*system).K,MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd((*system).K,MAT_FINAL_ASSEMBLY));

    // Remove constent kernel from pressure coefficient matrix
    Mat Me;
    MatCreate(PETSC_COMM_WORLD,&Me);
    ierr = MatSetSizes(Me,PETSC_DECIDE,PETSC_DECIDE,totalElem,totalElem);
    ierr = MatSetType(Me,MATMPIAIJ);
    ierr = MatSetUp(Me);

    for (int j=0; j<totalElem; j++){
    for (int i=0; i<totalElem; i++){
       MatSetValue(Me,j,i,1.0/(double)totalElem,INSERT_VALUES); 
    }}

    MatAssemblyBegin(Me,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(Me,MAT_FINAL_ASSEMBLY);

//    MatView(Me, PETSC_VIEWER_STDOUT_WORLD);

    MatAXPY(system->Cs,-1.0, Me,DIFFERENT_NONZERO_PATTERN);
    MatAXPY(system->Cd,-1.0, Me,DIFFERENT_NONZERO_PATTERN);

    PetscFunctionReturn(0);
}
