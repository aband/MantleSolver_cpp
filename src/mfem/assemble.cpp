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
    double theta = 0;
    double mu_s = 1.0;
    double mu_f = 1.0;
    double inv_k0 = 10e08;
    double rho_r = 2800/3300;
    double gx = 0;
    double gy = 10;

    double phi_f = 0.0;
    double phi_s = 1.0;

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
        for (unsigned int i=0; i<8; i++){(*locmatrix).ad[j*8+i] = 0.0;}}

    for (unsigned int j=0; j<12; j++){(*locmatrix).bs[j] = 0.0;
                                      (*locmatrix).rhs[j] = 0.0;
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

        phi_f = AssignPorosity(mapped, (*physproperty).l);
        phi_s = 1 - phi_f;

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

                (*locmatrix).as[i+j*12] += gw*jac* 2*mu_s*phi_s * (A1*A2+B1*B2*2+C1*C2);
            }

            (*locmatrix).bs[j] += gw*jac*div1 * 1;

            (*locmatrix).rhs[j] += gw*jac*(1-phi_f)*rho_r*(gx*brwork[j][0] + 
                                                           gy*brwork[j][1]);

        }

        std::array<vertex, 8> hdivwork = hdiv_.ComputeHdivmixed(basis_,mapped);

        for (unsigned int j=0; j<8; j++){
            for (unsigned int i=0; i<8; i++){
                (*locmatrix).ad[i+j*8] += gw*jac*mu_f*inv_k0* 
                                (hdivwork[j][0]*hdivwork[i][0] + 
                                 hdivwork[j][1]*hdivwork[i][1]);
            }
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
                                         Matrix * matrix,
                                         Vec * source){
 
    PetscErrorCode    ierr;
    PetscFunctionBeginUser;

    const valarray<double>& gwe = GaussWeightsEdge;
    const valarray<double>& gpe = GaussPointsEdge;
    const valarray<double>& gwf = GaussWeightsFace;
    const vector<vertex>& gpf = GaussPointsFace;

    PetscCall(MatCreate(PETSC_COMM_WORLD,&(*matrix).Ad));
    PetscCall(MatCreate(PETSC_COMM_WORLD,&(*matrix).As));
    PetscCall(MatCreate(PETSC_COMM_WORLD,&(*matrix).Bs));
    PetscCall(MatCreate(PETSC_COMM_WORLD,&(*matrix).Bd));
    PetscCall(MatCreate(PETSC_COMM_WORLD,&(*matrix).Cs));
    PetscCall(MatCreate(PETSC_COMM_WORLD,&(*matrix).Cd));
    PetscCall(MatCreate(PETSC_COMM_WORLD,&(*matrix).K));

    double *localrhs;

    PetscCall(VecGetArray(*source, &localrhs));

    hdiv_.ComputeTotalDOF(mi);

    // Set zeros to right hand side vector
    for (unsigned int k=0; k<br_.getDOF(); k++){localrhs[k] = 0.0;}

    int totalElem = mi.MPIglobalCellSize[0] * mi.MPIglobalCellSize[1];

    PetscCall(MatSetSizes((*matrix).Ad,PETSC_DECIDE,PETSC_DECIDE,
                                       hdiv_.getDOF(),hdiv_.getDOF()));
    PetscCall(MatSetSizes((*matrix).As,PETSC_DECIDE,PETSC_DECIDE,
                                       br_.getDOF(),br_.getDOF()));

    PetscCall(MatSetSizes((*matrix).Bd,PETSC_DECIDE,PETSC_DECIDE,
                                       hdiv_.getDOF(),totalElem));
    PetscCall(MatSetSizes((*matrix).Bs,PETSC_DECIDE,PETSC_DECIDE,
                                       br_.getDOF(),totalElem));

    PetscCall(MatSetSizes((*matrix).Cd,PETSC_DECIDE,PETSC_DECIDE,
                                       totalElem,totalElem));
    PetscCall(MatSetSizes((*matrix).Cs,PETSC_DECIDE,PETSC_DECIDE,
                                       totalElem,totalElem));

    PetscCall(MatSetSizes((*matrix).K,PETSC_DECIDE,PETSC_DECIDE,
                                       totalElem,totalElem));

    PetscCall(MatSetType((*matrix).Ad,MATMPIAIJ));
    PetscCall(MatSetType((*matrix).As,MATMPIAIJ));
    PetscCall(MatSetType((*matrix).Bd,MATMPIAIJ));
    PetscCall(MatSetType((*matrix).Bs,MATMPIAIJ));
    PetscCall(MatSetType((*matrix).Cd,MATMPIAIJ));
    PetscCall(MatSetType((*matrix).Cs,MATMPIAIJ));
    PetscCall(MatSetType((*matrix).K,MATMPIAIJ));
 
    ierr = MatSetUp((*matrix).As);CHKERRQ(ierr);
    ierr = MatSetUp((*matrix).Bs);CHKERRQ(ierr);
    ierr = MatSetUp((*matrix).Bd);CHKERRQ(ierr);
    ierr = MatSetUp((*matrix).Cs);CHKERRQ(ierr);
    ierr = MatSetUp((*matrix).Cd);CHKERRQ(ierr);
    ierr = MatSetUp((*matrix).K);CHKERRQ(ierr);
    ierr = MatSetUp((*matrix).Ad);CHKERRQ(ierr);

    LocMatrix * locmatrix = (LocMatrix *)malloc(sizeof(LocMatrix));

    int checkSizeM, checkSizeN;
    MatGetSize((*matrix).Ad, &checkSizeM, &checkSizeN);

    for (unsigned int n=0; n<totalElem; n++){
        indice globalElem = Bend(mi,n);

        basis_.GetCorners(mi,globalElem); 

        std::array<int, 8> tmp = hdiv_.LocalToGlobal(mi,globalElem);

        // Compute local matrix
        AssignLocMatrix(mi,basis_,hdiv_,br_,locmatrix,physproperty,gwe,gpe,gwf,gpf);

        const int ISDarcy[8] = 
           {tmp[0],tmp[1],tmp[2],tmp[3],tmp[4],tmp[5],tmp[6],tmp[7]};

        const double Bdv[8] = {(*locmatrix).bd[0],(*locmatrix).bd[1],
                               (*locmatrix).bd[2],(*locmatrix).bd[3],
                               (*locmatrix).bd[4],(*locmatrix).bd[5],
                               (*locmatrix).bd[6],(*locmatrix).bd[7]}; 

        const int idxm = n;

        PetscCall(MatSetValuesBlocked((*matrix).Bd, 8, ISDarcy, 1, &idxm, Bdv, ADD_VALUES));

        for (unsigned int l=0; l<8; l++){
            const double Adv[8] = {(*locmatrix).ad[0+l*8],(*locmatrix).ad[1+l*8],
                                   (*locmatrix).ad[2+l*8],(*locmatrix).ad[3+l*8],
                                   (*locmatrix).ad[4+l*8],(*locmatrix).ad[5+l*8],
                                   (*locmatrix).ad[6+l*8],(*locmatrix).ad[7+l*8]};
            const int idxm = ISDarcy[l];
            PetscCall(MatSetValuesBlocked((*matrix).Ad,1,&idxm,8,ISDarcy,Adv,ADD_VALUES));
        }

        PetscCall(MatSetValue((*matrix).Cd,n,n,(*locmatrix).cd,ADD_VALUES));

        // Assemble Stokes part

        const std::array<int, 12> tmp2 = br_.LocalToGlobal(mi,globalElem);

        const int ISStokes[12] = {tmp2[0],tmp2[1],tmp2[2],tmp2[3],
                                  tmp2[4],tmp2[5],tmp2[6],tmp2[7],
                                  tmp2[8],tmp2[9],tmp2[10],tmp2[11]};

        for (unsigned int k=0; k<12; k++){
            localrhs[ISStokes[k]] += (*locmatrix).rhs[k];
        }

        const double Bsv[12] = {(*locmatrix).bs[0],(*locmatrix).bs[1],
                                (*locmatrix).bs[2],(*locmatrix).bs[3],
                                (*locmatrix).bs[4],(*locmatrix).bs[5],
                                (*locmatrix).bs[6],(*locmatrix).bs[7],
                                (*locmatrix).bs[8],(*locmatrix).bs[9],
                                (*locmatrix).bs[10],(*locmatrix).bs[11]};

        const int idxms = n;
        PetscCall(MatSetValuesBlocked((*matrix).Bs,12,ISStokes,1,&idxms,Bsv,ADD_VALUES));

        for (unsigned int l=0; l<12; l++){
            const double Asv[12] = {(*locmatrix).as[0+l*12], (*locmatrix).as[1+l*12],
                                    (*locmatrix).as[2+l*12], (*locmatrix).as[3+l*12],
                                    (*locmatrix).as[4+l*12], (*locmatrix).as[5+l*12],
                                    (*locmatrix).as[6+l*12], (*locmatrix).as[7+l*12],
                                    (*locmatrix).as[8+l*12], (*locmatrix).as[9+l*12],
                                    (*locmatrix).as[10+l*12], (*locmatrix).as[11+l*12]};

            const int idxm = ISStokes[l];
            PetscCall(MatSetValuesBlocked((*matrix).As, 1, &idxm, 12, ISStokes, Asv, ADD_VALUES));
        }

        PetscCall(MatSetValue((*matrix).Cs, n, n, (*locmatrix).cs, ADD_VALUES));

        PetscCall(MatSetValue((*matrix).K, n, n, (*locmatrix).k, ADD_VALUES));

    }

    PetscCall(VecRestoreArray(*source,&localrhs));

    // Assemble all block matrix
    PetscCall(MatAssemblyBegin((*matrix).As,MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd((*matrix).As,MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyBegin((*matrix).Ad,MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd((*matrix).Ad,MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyBegin((*matrix).Bs,MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd((*matrix).Bs,MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyBegin((*matrix).Bd,MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd((*matrix).Bd,MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyBegin((*matrix).Cs,MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd((*matrix).Cs,MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyBegin((*matrix).Cd,MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd((*matrix).Cd,MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyBegin((*matrix).K,MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd((*matrix).K,MAT_FINAL_ASSEMBLY));

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

    MatAXPY(matrix->Cs,-1.0, Me,DIFFERENT_NONZERO_PATTERN);
    MatAXPY(matrix->Cd,-1.0, Me,DIFFERENT_NONZERO_PATTERN);

    return ierr;
}
