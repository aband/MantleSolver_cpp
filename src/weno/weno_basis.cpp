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

// This is the correct function to use
void WenoPrepare::CreateSmoothnessIndicator(const WenoMesh*& wm, double gamma, double Theta, double omega_l){

    omega_0 = omega_l;

    for (auto & cell: index_set_stencil){
        int target_i = target_cell[0]-wm->ghost;
        int target_j = target_cell[1]-wm->ghost;
        sigma2 += pow(wm->lsol[target_j][target_i] - wm->lsol[target_j+cell[1]][target_i+cell[0]],2);
    }

    sigma2 *= 1.0/(double)(index_set_stencil.size()-1);
    sigma2 = pow(sigma2,gamma/2.0);

    omega_hat = omega_l/(sigma2+epsilon_0*pow(h,Theta));

    Theta_l = Theta;
    gamma_l = gamma;
}

void WenoPrepare::CreateSmoothnessIndicator(const WenoMesh*& wm, int c, double gamma, double Theta, double omega_l){

    omega_0 = omega_l;

    index_set temp_index_set;

    temp_index_set.assign(index_set_stencil.begin(),index_set_stencil.end());

    // Erase center cell out of index set
    temp_index_set.erase(temp_index_set.begin()+c);

    int totali = wm->M + 2*wm->ghost; 

    for (auto & cell: temp_index_set){
        int target_i = target_cell[0]-wm->ghost;
        int target_j = target_cell[1]-wm->ghost;

        point p0 = {0.0,0.0};
        point p1 = {0.0,0.0};
        for (auto & c: corner_index){
            int t0 = target_cell[0]+c[0];
            int t1 = target_cell[1]+c[1];

            p0 += wm->lmesh[t1*totali+t0]/4.0;

            t0 = t0+cell[0];
            t1 = t1+cell[1];

            p1 += wm->lmesh[t1*totali+t0]/4.0;
        }

        p0 = (p0-p1)*(p0-p1);

        double work = wm->lsol[target_j][target_i] - wm->lsol[target_j+cell[1]][target_i+cell[0]];

        sigma2 += pow(work*h,2)/p0.sum();
    }
    sigma2 *= 1.0/(double)(index_set_stencil.size()-1);
    sigma2 = pow(sigma2,gamma);

    omega_hat = omega_l/(sigma2+epsilon_0*pow(h/10.0,Theta));

    Theta_l = Theta;
    gamma_l = gamma;
}

/*
 *Construct a object containing necessary information for a
 *single reconstruction at a target cell.
 */
WenoReconst::WenoReconst(point_index& target, const WenoMesh*& wm,
                         index_set& StencilLargeIndex, vector<int>& input_Lorder,
                         vector<index_set>& StencilSmallIndex, vector<int>& input_Sorder,
                         int center_indexl, vector<int> center_indexs){

    target_cell = target;

    Lorder = input_Lorder; 
    Sorder = input_Sorder; 

    r = max(Lorder[0],Lorder[1]);
    s = max(Sorder[0],Sorder[1]);

    StencilLarge = StencilLargeIndex;    
    StencilSmall = StencilSmallIndex;

    cil = center_indexl;
    cis = center_indexs;

    Theta = (double)s;
    gamma = Theta + (double)s;

    etas = ceil((max(r-s,s)+Theta)/2.0);
    etal = ceil((r-s+Theta)/2.0);

    lwp = new WenoPrepare(StencilLarge, target, Lorder);
    lwp->SetUpStencil(wm);
    lwp->CreateBasisCoeff(wm);
    //lwp->CreateSmoothnessIndicator(wm,etal,Theta);
    lwp->CreateSmoothnessIndicator(wm,r/2+1,r,0.5);
    //lwp->CreateSmoothnessIndicator(wm,cil,r/2+1,r,0.5);

    swp = new wpPtr[StencilSmall.size()]; 

    for (int e=0; e<StencilSmall.size(); e++){
        swp[e] = new WenoPrepare(StencilSmall[e], target, Sorder);
        swp[e]->SetUpStencil(wm);
        swp[e]->CreateBasisCoeff(wm);
        //swp[e]->CreateSmoothnessIndicator(wm,etas,Theta);
        swp[e]->CreateSmoothnessIndicator(wm,s/2+1,s,0.125);
        //swp[e]->CreateSmoothnessIndicator(wm,cis[e],s/2+1,s,0.125);
    }

}

WenoReconst::WenoReconst(point_index& target, const WenoMesh*& wm,
                         index_set& StencilLargeIndex, vector<int>& input_Lorder,
                         vector<index_set>& StencilSmallIndex, vector<int>& input_Sorder,
                         int center_indexl, vector<int> center_indexs,
                         WenoReconst*& wr){

    target_cell = target;

    Lorder = input_Lorder; 
    Sorder = input_Sorder; 

    r = max(Lorder[0],Lorder[1]);
    s = max(Sorder[0],Sorder[1]);

    StencilLarge = StencilLargeIndex;    
    StencilSmall = StencilSmallIndex;

    etas = ceil((max(r-s,s)+Theta)/2.0);
    etal = ceil((r-s+Theta)/2.0);

    cil = center_indexl;
    cis = center_indexs;

    Theta = (double)s;
    gamma = Theta + (double)s; 

    lwp = new WenoPrepare(StencilLarge, target, Lorder);
    lwp->SetUpStencil(wm);
    lwp->wenobasiscoeff = wr->lwp->wenobasiscoeff;
    //lwp->CreateSmoothnessIndicator(wm,etal,Theta);
    //lwp->CreateSmoothnessIndicator(wm,r/2+1,r,0.5);
    lwp->CreateSmoothnessIndicator(wm,cil,r/2+1,r,0.5);

    swp = new wpPtr[StencilSmall.size()]; 

    for (int e=0; e<StencilSmall.size(); e++){
        swp[e] = new WenoPrepare(StencilSmall[e], target, Sorder);
        swp[e]->SetUpStencil(wm);
        swp[e]->wenobasiscoeff = wr->swp[e]->wenobasiscoeff;
        //swp[e]->CreateSmoothnessIndicator(wm,etas,Theta);
        //swp[e]->CreateSmoothnessIndicator(wm,s/2+1,s,0.125);
        swp[e]->CreateSmoothnessIndicator(wm,cis[e],s/2+1,s,0.125);
    }

}

/*
 *void WenoReconst::WenoUpdate(const WenoMesh*& wm){
 *
 *    //lwp->CreateSmoothnessIndicator(wm,etal,Theta);
 *    lwp->CreateSmoothnessIndicator(wm,gamma,Theta,0.5);
 *
 *    for (int s=0; s<StencilSmall.size(); s++){
 *        //swp[s]->CreateSmoothnessIndicator(wm,etas,Theta);
 *        swp[s]->CreateSmoothnessIndicator(wm,gamma,Theta,0.125);
 *    }
 *
 *}
 *
 */
void WenoReconst::CreateWeights(){

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

void WenoReconst::CreateNewWeights2(){

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
    for (int k=0; k<StencilSmall.size(); k++){
        sweight[k] = omega_tilde[k+1];
    }
    lweight = omega_tilde[0];
}

void WenoReconst::CheckBasisCoeff(){
    printf("Basis polynomial on large stencil : \n");
    lwp->PrintBasisCoeff();
    printf("Basis polynomial on Small stencil : \n");
    for(int s=0; s<StencilSmall.size(); s++){
        swp[s]->PrintBasisCoeff();
        cout << endl;
    }
} 

void WenoReconst::CheckWeights(){
    printf("Check reconstruction weights : \n");
    printf("Large stencil weight: %f \n",lweight); 
    for (int i=0; i<StencilSmall.size(); i++){
        printf("Small stencil weights %f \n",sweight[i]);
    }
}

solution WenoReconstStencil(vector<int>& order, point_index& target, point target_point,
                            WenoPrepare*& wp,const WenoMesh*& wm){

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
