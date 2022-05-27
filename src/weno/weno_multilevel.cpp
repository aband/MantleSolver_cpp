#include "../../include/weno_multilevel.h"

/*
 *A more general way of handling multilevel weno reconstruction.
 */

double constfunc(valarray<double>& point,const vector<double>& param){
    return 1.0;
}

double poly(valarray<double>& point, const vector<int>& param){
    return pow(point[0],param[0])*pow(point[1],param[1]);
}

/*
 *Auxiliary functions during the calculation process.
 */

int factorial(int top, int bottom){
    assert(top>bottom || top==bottom);
    if (top==bottom){ return 1;}
    else{ return top*factorial(top-1,bottom);}
}

int factorial(int top){
    assert(top>0 || top==0);
    if (top==0){ return 1;}
    else{ return top*factorial(top-1);}
}

int * polynDeriv(int degree, int high){

    assert(degree < high || degree == high);

    int * multiplier = new int[high]();

    for (int i=0; i<high-degree; i++){
        multiplier[i] = factorial(degree+i,i);
    }

    return multiplier;
}

int * polynDeriv(int xdegree, int ydegree, int xhigh, int yhigh){

    int * multiplier = new int[xhigh*yhigh]();

    int * dxPloyn = polynDeriv(xdegree,xhigh);
    int * dyPloyn = polynDeriv(ydegree,yhigh);

    for (int j=0; j<yhigh; j++){
    for (int i=0; i<xhigh; i++){
        multiplier[i+j*yhigh] = dxPloyn[i]*dyPloyn[j];
    }}

    return multiplier;
}

// =========================================================================================
WenoStencil::WenoStencil(MeshInfo* mi,const int rangex[2], point_index& target){

    assert(target.size() == 1);
	 // Shift target point index with ghost region
    int targetx = target[0] -1 + mi->ghost_vertx[0];

    int orderx = rangex[1]-rangex[0]+1;

    polyn_order.push_back(orderx);

    assert(polyn_order.size() == target.size());

    stencil_size = orderx;

    myStencil <point_index> stencil(rangex[1]-rangex[0]+1);

    stencil.CreateStencil();

    for (int i=0; i<orderx; i++){
        valarray<int> start {rangex[0]};
        valarray<int> incre {i};
        valarray<int> target_shift {targetx};
        stencil(i) = start + incre + target_shift; 
    }

    stencil_index_set = stencil.GetStencil(); 

    // =======================================================
    stencil_center = {0.0};

    int fullx = mi->localsize[0]+2*mi->ghost_vertx[0];

    for (auto & c: mi->corner_index_1D){
        int flat = targetx + c[0];
        point corner = mi->lmesh[flat];
        center_cell_corners.push_back(corner);
        stencil_center += corner/(double)mi->corner_index_1D.size();
    }

    //h = pow(NumIntegralFace(center_cell_corners,{0.0},{0.0,0.0},1.0,constfunc),0.5);
}

WenoStencil::WenoStencil(MeshInfo* mi,const int rangex[2], const int rangey[2], point_index& target){

    assert(target.size() == 2);
	 // Shift target point index with ghost region
    int targetx = target[0] -1 + mi->ghost_vertx[0];
    int targety = target[1] -1 + mi->ghost_vertx[1];

    int orderx = rangex[1]-rangex[0]+1;
    int ordery = rangey[1]-rangey[0]+1;

    polyn_order.push_back(orderx);
    polyn_order.push_back(ordery);

    assert(polyn_order.size() == target.size());

    stencil_size = orderx*ordery;

    myStencil <point_index> stencil(orderx,ordery);

    stencil.CreateStencil();

    for (int j=0; j<ordery; j++){
    for (int i=0; i<orderx; i++){
        valarray<int> start {rangex[0],rangey[0]};
        valarray<int> incre {i,j};
        valarray<int> target_shift {targetx,targety};
        stencil(i,j) = start + incre + target_shift; 
    }}

    stencil_index_set = stencil.GetStencil(); 

    // =======================================================
    stencil_center = {0.0,0.0};

    int fullx = mi->localsize[0]+2*mi->ghost_vertx[0];

    for (auto & c: mi->corner_index_2D){
        int flat = (targety+c[1])*fullx+targetx+c[0];
        point corner = mi->lmesh[flat];
        center_cell_corners.push_back(corner);
        stencil_center += corner/(double)mi->corner_index_2D.size();
    }

    h = pow(NumIntegralFace(center_cell_corners,{0.0},{0.0,0.0},1.0,constfunc),0.5);

    CreateBasisPolyn(mi); 
    CreateSigma(mi);;
}

WenoStencil::WenoStencil(MeshInfo* mi,const int rangex[2], const int rangey[2], const int rangez[2], point_index& target){

    assert(target.size() == 3);
	 // Shift target point index with ghost region
    int targetx = target[0] -1 + mi->ghost_vertx[0];
    int targety = target[1] -1 + mi->ghost_vertx[1];
    int targetz = target[2] -1 + mi->ghost_vertx[2];

    int orderx = rangex[1]-rangex[0]+1;
    int ordery = rangey[1]-rangey[0]+1;
    int orderz = rangez[1]-rangez[0]+1;

    polyn_order.push_back(orderx);
    polyn_order.push_back(ordery);
    polyn_order.push_back(orderz);

    assert(polyn_order.size() == target.size());

    stencil_size = orderx*ordery*orderz;

    myStencil <point_index> stencil(orderx,ordery,orderz);

    stencil.CreateStencil();

    for (int k=0; k<orderz; k++){
    for (int j=0; j<ordery; j++){
    for (int i=0; i<orderx; i++){
        valarray<int> start {rangex[0],rangey[0],rangez[0]};
        valarray<int> incre {i,j,k};
        valarray<int> target_shift {targetx,targety,targetz};
        stencil(i,j,k) = start + incre + target_shift; 
    }}}

    stencil_index_set = stencil.GetStencil(); 

    // =======================================================
    stencil_center = {0.0,0.0,0.0};

    int fullx = mi->localsize[0]+2*mi->ghost_vertx[0];
    int fully = mi->localsize[1]+2*mi->ghost_vertx[1];

    for (auto & c: mi->corner_index_3D){
        int flat = (targetz+c[2])*fully*fullx+(targety+c[1])*fullx+targetx+c[0];
        point corner = mi->lmesh[flat];
        center_cell_corners.push_back(corner);
        stencil_center += corner/(double)mi->corner_index_3D.size();
    }

    //h = pow(NumIntegralFace(center_cell_corners,{0.0},{0.0,0.0},1.0,constfunc),0.5);

}

void WenoStencil::CheckWenoStencil(){

    for (int k=0; k<stencil_size; k++){
        cout << "(" ;
        for (int p=0; p<stencil_index_set[k].size(); p++){
            printf("%d ",stencil_index_set[k][p]);
        }
        cout << ") " ;
    }printf("\n");
}

// 1D Case ==========================================================

// ==================================================================

// 2D Case ==========================================================
void WenoStencil::CreateBasisPolyn(MeshInfo* mi){

    int fullx = mi->localsize[0]+2*mi->ghost_vertx[0];

    //Set up linear system for solving Basis coefficients.
    lapack_int n    = polyn_order.at(0)*polyn_order.at(1);
    lapack_int nrhs = n;
    lapack_int lda  = n;
    lapack_int ldb  = nrhs;

    double * a = new double [n*n];
    double * b = new double [n*nrhs];
    lapack_int * p = new int [n];

    for (int cell=0; cell<stencil_index_set.size(); cell++){
        for (int ypow=0; ypow<polyn_order.at(1); ypow++){
        for (int xpow=0; xpow<polyn_order.at(0); xpow++){
            vector< point > work;
            for (auto & c: mi->corner_index_2D){
                int stencilj = stencil_index_set[cell][1] + c[1];
                int stencili = stencil_index_set[cell][0] + c[0];
                work.push_back(mi->lmesh[stencilj*fullx+stencili]);
            }
            a[cell*n+ypow*polyn_order[0]+xpow] = NumIntegralFace(work,{xpow,ypow},stencil_center,h,poly);
        }}
    }
    fill(b,b+n*nrhs,0);
    for (int i=0; i<nrhs; i++){b[i*n+i]=a[n*i];}

    int err = LAPACKE_dgesv(LAPACK_ROW_MAJOR, n, nrhs, a, lda, p, b, ldb);
    if (err){
        printf("ERROR: Weno Basis Coefficient for order %d, %d. Error type %d \n",
                 polyn_order[0],polyn_order[1],err);
    }

    polyn = b;

}

double * WenoStencil::CreateBasisPolynDeriv(int xdegree, int ydegree, int k){

    assert(k<stencil_size);

    int * polynderiv = polynDeriv(xdegree,ydegree,polyn_order[0],polyn_order[1]); 

    double * product = new double [stencil_size]();

    for (int j=ydegree; j<polyn_order[1]; j++){
    for (int i=xdegree; i<polyn_order[0]; i++){
        int shiftj = j-ydegree;
        int shifti = i-xdegree;
        product[shiftj*(polyn_order[0])+shifti] = polyn[k*stencil_size + i+j*polyn_order[0]];
    }}

    for (int j=0; j<polyn_order[1]; j++){
    for (int i=0; i<polyn_order[0]; i++){
        product[i+j*polyn_order[0]] *= polynderiv[i+j*polyn_order[0]]; 
    }}

    delete polynderiv;

    return product;
}

void WenoStencil::CreateSigma(MeshInfo*& mi){

    CreateBasisPolynDeriv(1,0,0);

    double * sum = new double [stencil_size*stencil_size]();



/*
 *    for (int q = 0; q<polyn_order[1]; q++){
 *    for (int p = 0; p<polyn_order[0]; p++){
 *        double * Dpq = CreateBasisPolynDeriv(p,q,p,q);
 *        for (int q_prime=0; q_prime<polyn_order[1]; q_prime++){
 *        for (int p_prime=0; p_prime<polyn_order[0]; p_prime++){
 *            double * Dpq_prime = CreateBasisPolynDeriv(p_prime,q_prime,p_prime,q_prime);
 *        }}
 *    }}
 *
 */
}


// 3D Case ==========================================================

// ==================================================================

void WenoStencil::PrintBasisPolyn(){
    int n = polyn_order[0]*polyn_order[1];
    for (int j=0; j<n; j++){
    for (int i=0; i<n; i++){
        printf("coeff = %.4f ",polyn[i+j*n]);
    }cout<<endl;}
}
