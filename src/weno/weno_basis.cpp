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
        printf("(%.2f,%.2f)  ",c[0],c[1]);
    } cout<<endl;
    printf("The center of the target cell is (%.2f,%.2f) \n",center[0],center[1]);

    printf("The area is %.3f \n",h*h);

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
        sigma += pow(wm->lsol[target_j][target_i] - wm->lsol[target_j+cell[1]][target_i+cell[0]],2) ;
    }

    sigma *= 1.0/(double)(index_set_stencil.size()-1);

    sigma = pow(sigma, eta);

    omega = 1.0/pow(sigma+epsilon_0*h,(double)max(polynomial_order[0],polynomial_order[1])/Theta);
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

    int target_i = target[0] - wm->ghost;
    int target_j = target[1] - wm->ghost;

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

    double sum = accumulate(omega.begin(), omega.end(), decltype(omega)::value_type(0));
    sum = sum/omega_m;

    vector <double> sweight;
    for (int i=0; i<swp.size(); i++){
        sweight.push_back(omega[i]/omega_m*(1-lwp->omega/omega_m)/sum);
    }

    double lweight = 1.0 - accumulate(sweight.begin(), sweight.end(), decltype(sweight)::value_type(0));

    // Now calculate reconstruction with caluclated weightings.
    reconst = 0.0; 

    /*
     *Calcualte basis function evaluated at the given target point.
     */
    int sn = Sorder[0]*Sorder[1];
    for (auto & wp: swp){

        double work = 0.0;
        for (int p=0; p<sn; p++){
            for(int ypow=0; ypow<Sorder[1]; ypow++){
            for(int xpow=0; xpow<Sorder[0]; xpow++){
                int o = ypow*Sorder[0]+xpow;
                work += wm->lsol[target_j][target_i]*
                        wp->wenobasiscoeff[o*sn+p]*
                        poly(target_point,{xpow,ypow});
            }}
        }
        swork.push_back(work);

    }

    int ln = Lorder[0]*Lorder[1];
    for (int p=0; p<ln; p++){
        for(int ypow=0; ypow<Lorder[1]; ypow++){
        for(int xpow=0; xpow<Lorder[0]; xpow++){
            int o = ypow*Lorder[0]+xpow;
            lwork += wm->lsol[target_j][target_i]*
                     lwp->wenobasiscoeff[o*sn+p]*
                     poly(target_point,{xpow,ypow});
        }}
    }

    // Calculate the reconstruction value
    for (int i=0; i<sweight.size(); i++){
        reconst += sweight.at(i)*swork.at(i); 
    }
    reconst += lweight*lwork;

    return reconst;
}
