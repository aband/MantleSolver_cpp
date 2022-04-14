#include "../../include/weno.h"

double constfunc(valarray<double>& point,const vector<double>& param){
    return 1.0;
}

double poly(valarray<double>& point, const vector<int>& param){
    return pow(point[0],param[0])*pow(point[1],param[1]);
}

/*
 *The following member functions are invariant to different index, and
 *reconstruction order.
 */
void WenoStencil::SetUpStencil(const WenoMesh*& wm){

    // Unpack parameters carried by WenoMesh
    int totali = wm->M+2*wm->ghost;

    center = {0.0,0.0};
    /*
     *Locate target cell four corners.
     *Calculate center point at the same time.
     */
    for (auto & c: corner_index){
        point corner =  wm->lmesh[(target_cell[1]+c[1])*totali+target_cell[0]+c[0]]; 
        target_cell_corners.push_back(corner);
        center += corner/4.0;
    }

    h = pow(NumIntegralFace(target_cell_corners,{0.0},{0.0,0.0},1.0,constfunc),0.5);
}

void WenoStencil::PrintSingleStencil(){

    printf("The target cell index is (%d,%d) \n", target_cell[0], target_cell[1]);
    printf("The four corners are : \n");
    for (auto & c: target_cell_corners){
        printf("(%.5f,%.5f)  ",c[0],c[1]);
    } cout<<endl;
    printf("The center of the target cell is (%.5f,%.5f) \n",center[0],center[1]);

    printf("The area is %.12f \n",h*h);

}

void WenoPrepare::CreateBasisCoeff(const WenoMesh*& wm){
  
     int totali = wm->M+2*wm->ghost;

    /*
     *Set up linear system for solving Basis coefficients.
     */
     lapack_int n    = polynomial_order[0]*polynomial_order[1];
     lapack_int nrhs = n;
     lapack_int lda  = n;
     lapack_int ldb  = nrhs;

     double * a = new double [n*n];
     double * b = new double [n*nrhs];
     lapack_int * p = new int [n];
     
     for (int cell=0; cell<index_set_stencil.size(); cell++){
         for (int ypow=0; ypow<polynomial_order[1]; ypow++){
         for (int xpow=0; xpow<polynomial_order[0]; xpow++){
             cell_corners work;
             for (auto & c: corner_index){
                 int stencilj = target_cell[1]+index_set_stencil[cell][1]+c[1];
                 int stencili = target_cell[0]+index_set_stencil[cell][0]+c[0];

                 work.push_back(wm->lmesh[stencilj*totali+stencili]);
             }
             a[cell*n+ypow*polynomial_order[0]+xpow] = NumIntegralFace(work,{xpow,ypow},center,h,poly);
         }}
     }
     fill(b,b+n*nrhs,0);
     for (int i=0; i<nrhs; i++){b[i*n+i]=a[n*i];}

     int err = LAPACKE_dgesv(LAPACK_ROW_MAJOR, n, nrhs, a, lda, p, b, ldb);
     if (err){
         printf("ERROR: Weno Basis Coefficient for order %d, %d. Error type %d \n",
                 polynomial_order[0],polynomial_order[1],err);
     }

     wenobasiscoeff = b;
}

void WenoPrepare::PrintBasisCoeff(){
    int n = polynomial_order[0]*polynomial_order[1];
    for (int j=0; j<n; j++){
    for (int i=0; i<n; i++){
        printf("coeff = %.4f ",wenobasiscoeff[i*n+j]);
    }cout<<endl;}
}

void WenoPrepare::CreateSmoothnessIndicator(const WenoMesh*& wm, double eta, double Theta){

    for (auto & cell: index_set_stencil){
        int target_i = target_cell[0]-wm->ghost;
        int target_j = target_cell[1]-wm->ghost;
        sigma += pow(wm->lsol[target_j][target_i] - wm->lsol[target_j+cell[1]][target_i+cell[0]],2);
    }

    sigma *= 1.0/(double)(index_set_stencil.size()-1);

    sigma = pow(sigma, eta);

    omega = 1.0/pow(sigma+epsilon_0*pow(h,Theta),(double)max(polynomial_order[0],polynomial_order[1])/Theta);
}

void WenoPrepare::CreateNewSmoothnessIndicator(const WenoMesh*& wm, double gamma, double Theta){

    for (auto & cell: index_set_stencil){
        int target_i = target_cell[0]-wm->ghost;
        int target_j = target_cell[1]-wm->ghost;
        sigma += pow(wm->lsol[target_j][target_i] - wm->lsol[target_j+cell[1]][target_i+cell[0]],2);
    }
    sigma *= 1.0/(double)(index_set_stencil.size()-1);
    sigma = pow(sigma,gamma/2.0);

    omega_hat = omega_0/(sigma+epsilon_0*pow(h,Theta));
}

/*
 *Construct a object containing necessary information for a
 *single reconstruction at a target cell.
 */
WenoReconst::WenoReconst(point_index& target, const WenoMesh*& wm,
                         index_set& StencilLargeIndex, vector<int>& input_Lorder,
                         vector<index_set>& StencilSmallIndex, vector<int>& input_Sorder){

    // Restore parameters inside the class
    target_cell = target;

    Lorder = input_Lorder; 
    Sorder = input_Sorder; 

    StencilLarge = StencilLargeIndex;
    StencilSmall = StencilSmallIndex;

    r = max(Lorder[0],Lorder[1]);
    s = max(Sorder[0],Sorder[1]);

    etas = ceil((max(r-s,s)+Theta)/2.0);
    etal = ceil((r-s+Theta)/2.0);

    lwp = new WenoPrepare(StencilLarge, target, Lorder);
    lwp->SetUpStencil(wm);
    lwp->CreateBasisCoeff(wm);
    lwp->omega_0 = 0.5;
    //lwp->CreateSmoothnessIndicator(wm,etal,Theta);
    lwp->CreateNewSmoothnessIndicator(wm,8,Theta);

    swp = new wpPtr[StencilSmall.size()]; 

    for (int s=0; s<StencilSmall.size(); s++){
        swp[s] = new WenoPrepare(StencilSmall[s], target, Sorder);
        swp[s]->SetUpStencil(wm);
        swp[s]->CreateBasisCoeff(wm);
        swp[s]->omega_0 = 0.125;
        //swp[s]->CreateSmoothnessIndicator(wm,etas,Theta);
        swp[s]->CreateNewSmoothnessIndicator(wm,8,Theta);
    }

}

void WenoReconst::CreateCoefficients(){

    vector<double> omega;
    omega.push_back(lwp->omega);
    for (int s=0; s<StencilSmall.size(); s++){
        omega.push_back(swp[s]->omega);
    }

    /*
     *Create weighting based on smoothness indicators
     */
    double omega_m = *max_element(omega.begin(),omega.end());    
    omega_m = max(omega_m,lwp->omega);

    double sum = accumulate(omega.begin(), omega.end(), decltype(omega)::value_type(0));
    sum = sum/omega_m;

    double sum_sweight = 0.0;

    sweight = new double[StencilSmall.size()];
    for (int k=0; k<StencilSmall.size(); k++){
        sweight[k] = omega[k+1]/omega_m*(1-lwp->omega/omega_m)/sum;
        sum_sweight += sweight[k];
    }

    lweight = 1.0 - sum_sweight;
}

void WenoReconst::CreateNewWeights(){
    vector<double> omega_tilde;

    omega_tilde.push_back(lwp->omega_hat);
    for (int s=0; s<StencilSmall.size(); s++){
        omega_tilde.push_back(swp[s]->omega_hat);
    }

    double sum = accumulate(omega_tilde.begin(), omega_tilde.end(), decltype(omega_tilde)::value_type(0));
    for (auto & o: omega_tilde){
        o /= sum;
    }
  
    sweight = new double[StencilSmall.size()];
    lweight = 1.0;
    for (int k=0; k<StencilSmall.size(); k++){
        sweight[k] = omega_tilde[k+1]*pow(1.0 - omega_tilde[0]/lwp->omega_0,2);
        lweight -= sweight[k];
    }
}

void WenoReconst::CheckBasisCoeff(){
    printf("Basis polynomial on large stencil : \n");
    lwp->PrintBasisCoeff();
} 

void WenoReconst::CheckWeights(){
    printf("The computed linear weights are: \n");
    printf("a0 = %.16f , a1 = %.16f , a2 = %.16f , a3 = %.16f , a4 = %.16f \n",lweight,sweight[0],sweight[1],sweight[2],sweight[3]);
}

solution WenoReconstStencil(vector<int>& order, point_index& target, point target_point,
                            WenoPrepare*& wp, const WenoMesh*& wm){

    solution reconst = 0.0;

    int n = order[0]*order[1];
    for (int p=0; p<n; p++){
        int target_i = target[0] - wm->ghost + wp->index_set_stencil[p][0];
        int target_j = target[1] - wm->ghost + wp->index_set_stencil[p][1];

        for (int ypow = 0; ypow<order[1]; ypow++){
        for (int xpow = 0; xpow<order[0]; xpow++){
            int o = ypow*order[0] + xpow;
            point target = (target_point - wp->center)/wp->h;
            reconst += wm->lsol[target_j][target_i] *
                       wp->wenobasiscoeff[o*n+p] *
                       poly(target,{xpow,ypow});
        }}
    }

    return reconst;
}

solution WenoReconst::PointReconstruction(const WenoMesh*& wm, point target_point){

    solution work = 0.0;

    vector <double> swork;

    for (int s=0; s<StencilSmall.size(); s++){
        swork.push_back(WenoReconstStencil(Sorder, target_cell, target_point, swp[s], wm));
    }

    double lwork = WenoReconstStencil(Lorder, target_cell, target_point, lwp, wm);

    for (int s=0; s<StencilSmall.size(); s++){
        work += sweight[s]*swork[s];
    }
    work += lweight*lwork;

    return work;
}

// Define a reconstruction method
solution WenoPointReconst(index_set& StencilLarge, vector<index_set>& StencilSmall,const WenoMesh*& wm,
                          point_index& target, vector<int>& Sorder, vector<int>& Lorder,
                          point target_point){

    solution reconst;

    /*
     *Create parameters and corresponding containers
     */
    double Theta = 2.0;
    
    vector< WenoPrepare* > swp;
    vector<double> omega;

    int r = max(Lorder[0],Lorder[1]);
    int s = max(Sorder[0],Sorder[1]);

    double etas = ceil((max(r-s,s)+Theta)/2.0);
    double etal = ceil((r-s+Theta)/2.0);

    vector <double> swork;
    double lwork;

    /*
     *For small stencil
     */
    for (auto & is: StencilSmall){
        WenoPrepare * wp = new WenoPrepare(is,target,Sorder);
        wp->SetUpStencil(wm);
        wp->CreateBasisCoeff(wm);

        wp->CreateSmoothnessIndicator(wm,etas,Theta);
        omega.push_back(wp->omega);
        swp.push_back(wp);
    }

    /*
     *For large stencils
     */
    WenoPrepare * lwp = new WenoPrepare(StencilLarge,target,Lorder);
    lwp->SetUpStencil(wm);
    lwp->CreateBasisCoeff(wm);
    lwp->CreateSmoothnessIndicator(wm,etal,Theta);

    /*
     *Create weighting based on smoothness indicators
     */
    double omega_m = *max_element(omega.begin(),omega.end());    
    omega_m = max(omega_m,lwp->omega);

    double sum = accumulate(omega.begin(), omega.end(), decltype(omega)::value_type(0));
    sum = sum/omega_m;

    vector <double> sweight;
    for (int i=0; i<swp.size(); i++){
        sweight.push_back(omega[i]/omega_m*(1-lwp->omega/omega_m)/sum);
    }

    double lweight = 1.0 - accumulate(sweight.begin(), sweight.end(), decltype(sweight)::value_type(0));

    // Linear weights
    //sweight = {0.25,0.25,0.25,0.25};
    //lweight = 0.0;

    // Now calculate reconstruction with caluclated weightings.
    reconst = 0.0; 

    /*
     *Calcualte basis function evaluated at the given target point.
     */
    for (auto & wp: swp){
        swork.push_back(WenoReconstStencil(Sorder,target,target_point,wp,wm));
    }

    lwork = WenoReconstStencil(Lorder,target,target_point,lwp,wm);

    // Calculate the reconstruction value
    for (int i=0; i<sweight.size(); i++){
        reconst += sweight.at(i)*swork.at(i); 
    }
    reconst += lweight*lwork;

    return reconst;
}

void WenoReconst::CreateQuasiDerivative(vector<int>& order, point_index& target, point target_point,
                                        int currentStencilIndex, 
                                        WenoPrepare*& wp, const WenoMesh*& wm, double weight){

    // This function creates quasi derivatives of the reconstruction values.
    // Setup offset i and j based on largeStencil index set
    int offseti = StencilLarge[0][0];
    int offsetj = StencilLarge[0][1]; 

    int N = Lorder[0];

    int n = order[0]*order[1];
    double phi = 0.0;

    for (int p=0; p<n; p++){
        int locali = wp->index_set_stencil[p][0] - offseti;
        int localj = wp->index_set_stencil[p][1] - offsetj; 
        
        int target_i = target[0]-wm->ghost + wp->index_set_stencil[p][0];
        int target_j = target[1]-wm->ghost + wp->index_set_stencil[p][1];

        for (int ypow = 0; ypow<order[1]; ypow++){
        for (int xpow = 0; xpow<order[0]; xpow++){
            int o = ypow*order[0] + xpow;
            point target = (target_point - wp->center)/wp->h;
            phi += wp->wenobasiscoeff[o*n+p]*poly(target,{xpow,ypow}) ;
        }}

        deriv[localj*N+locali] += weight*phi;
    }

}

/*
 *void WenoReconst::PrepareWenoDerivative(vector<int>){
 *
 *}
 */
