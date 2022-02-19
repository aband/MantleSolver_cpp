#include "../../include/newweno.h"

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

void WenoPrepare::CreateSmoothnessIndicator(const WenoMesh*& wm, double eta){

    for (auto & cell: index_set_stencil){
        int target_i = target_cell[0]-wm->ghost;
        int target_j = target_cell[1]-wm->ghost;
        sigma += pow(wm->lsol[target_j][target_i] - wm->lsol[target_j+cell[1]][target_i+cell[0]],2) ;
    }

    sigma *= 1.0/(double)(index_set_stencil.size()-1);

    sigma = pow(sigma, eta);
}

// Define a reconstruction method
solution WenoReconstruction(index_set& StencilLarge, vector<index_set>& StencilSmall, WenoMesh*& wm,
                            point_index target, const vector<int>& Sorder, const vector<int>& Lorder,
                            point target_point){

    solution reconst;

    // For small stencils
    vector< WenoPrepare* > SmallBasisCoeff;
    vector< WenoPrepare* > SmallSigma;


    return reconst;
}
