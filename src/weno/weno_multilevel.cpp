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
WenoStencil::WenoStencil(const MeshInfo& mi,const int rangex[2], point_index& target){

    assert(target.size() == 1);
	 // Shift target point index with ghost region
    int targetx = target[0] + mi.ghost_vertx[0];

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

    int fullx = mi.localsize[0]+2*mi.ghost_vertx[0];

    for (auto & c: mi.corner_index_1D){
        int flat = targetx + c[0];
        point corner = mi.lmesh[flat];
        center_cell_corners.push_back(corner);
        stencil_center += corner/(double)mi.corner_index_1D.size();
    }

    //h = pow(NumIntegralFace(center_cell_corners,{0.0},{0.0,0.0},1.0,constfunc),0.5);

    // Restore target cell and vertex index
    targetCell_ = target;
    targetVertx_ = {targetx};
}

WenoStencil::WenoStencil(const MeshInfo& mi,const int rangex[2], const int rangey[2], point_index& target){

    assert(target.size() == 2);
	 // Shift target point index with ghost region
    int targetx = target[0] + mi.ghost_vertx[0];
    int targety = target[1] + mi.ghost_vertx[1];

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

    int fullx = mi.localsize[0]+2*mi.ghost_vertx[0];

    for (auto & c: mi.corner_index_2D){
        int flat = (targety+c[1])*fullx+targetx+c[0];
        point corner = mi.lmesh[flat];
        center_cell_corners.push_back(corner);
        stencil_center += corner/(double)mi.corner_index_2D.size();
    }

    h = pow(NumIntegralFace(center_cell_corners,{0.0},{0.0,0.0},1.0,constfunc),0.5);

    // Restore target cell and vertex index
    targetCell_ = target;
    targetVertx_ = {targetx,targety};

    // new all pointers in constructor
    polyn = new double [stencil_size*stencil_size]();  
    sigma = new double [stencil_size*stencil_size]();
    polynderiv_ = new int * [stencil_size];

    for (int i=0; i<stencil_size; i++){
        polynderiv_[i] = new int [stencil_size]();
    }

    CreatePolynDerivMulti();
    CreateBasisPolyn(mi); 
    CreateSigma(mi);

}

WenoStencil::WenoStencil(const MeshInfo& mi,const int rangex[2], const int rangey[2], const int rangez[2], point_index& target){

    assert(target.size() == 3);
	 // Shift target point index with ghost region
    int targetx = target[0] + mi.ghost_vertx[0];
    int targety = target[1] + mi.ghost_vertx[1];
    int targetz = target[2] + mi.ghost_vertx[2];

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

    int fullx = mi.localsize[0]+2*mi.ghost_vertx[0];
    int fully = mi.localsize[1]+2*mi.ghost_vertx[1];

    for (auto & c: mi.corner_index_3D){
        int flat = (targetz+c[2])*fully*fullx+(targety+c[1])*fullx+targetx+c[0];
        point corner = mi.lmesh[flat];
        center_cell_corners.push_back(corner);
        stencil_center += corner/(double)mi.corner_index_3D.size();
    }

    //h = pow(NumIntegralFace(center_cell_corners,{0.0},{0.0,0.0},1.0,constfunc),0.5);

    // Restore target cell and vertex index
    targetCell_ = target;
    targetVertx_ = {targetx, targety, targetz};
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
void WenoStencil::CreateBasisPolyn(const MeshInfo& mi){

    int fullx = mi.localsize[0]+2*mi.ghost_vertx[0];

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
            for (auto & c: mi.corner_index_2D){
                int stencilj = stencil_index_set[cell][1] + c[1];
                int stencili = stencil_index_set[cell][0] + c[0];
                work.push_back(mi.lmesh[stencilj*fullx+stencili]);
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

//    polyn = b;

    // Deep copy b to polyn
    for (int i=0; i<n*nrhs; i++){
        polyn[i] = b[i];
    }


    delete [] a;
    delete [] b;
    delete [] p;
}

void WenoStencil::CreatePolynDerivMulti(){

    for (int ydegree=0; ydegree<polyn_order[1]; ydegree++){
    for (int xdegree=0; xdegree<polyn_order[0]; xdegree++){
        polynderiv_[ydegree*polyn_order[0]+xdegree] = polynDeriv(xdegree,ydegree,polyn_order[0],polyn_order[1]);
    }}

}

void WenoStencil::CheckPolynDerivMulti(){
     for (int ydegree=0; ydegree<polyn_order[1]; ydegree++){
     for (int xdegree=0; xdegree<polyn_order[0]; xdegree++){
         for (int yd=0; yd<polyn_order[1]; yd++){
         for (int xd=0; xd<polyn_order[0]; xd++){
             cout << polynderiv_[ydegree*polyn_order[0]+xdegree][yd*polyn_order[0]+xd] << "  ";
         }cout << endl;} cout << endl;
     }cout << endl;} cout << endl;
}

double * WenoStencil::CreateBasisPolynDeriv(int xdegree, int ydegree, int k){

    assert(k<stencil_size);

    //int * polynderiv = polynDeriv(xdegree,ydegree,polyn_order[0],polyn_order[1]); 

    int * polynderiv = polynderiv_[ydegree*polyn_order[0]+xdegree];

    double * product = new double [stencil_size]();

    for (int j=ydegree; j<polyn_order[1]; j++){
    for (int i=xdegree; i<polyn_order[0]; i++){
        int shiftj = j-ydegree;
        int shifti = i-xdegree;
        product[shiftj*(polyn_order[0])+shifti] = polyn[k + (i+j*polyn_order[0])*stencil_size];
    }}

    for (int j=0; j<polyn_order[1]; j++){
    for (int i=0; i<polyn_order[0]; i++){
        product[i+j*polyn_order[0]] *= polynderiv[i+j*polyn_order[0]]; 
    }}

    return product;
}

void WenoStencil::CreateSigma(const MeshInfo& mi){

    int fullx = mi.localsize[0]+2*mi.ghost_vertx[0];

    // Get integral domain
    vector< point > work; 
    for (auto & corner: mi.corner_index_2D){
        int vertxj = targetVertx_[1] + corner[1];
        int vertxi = targetVertx_[0] + corner[0];
        work.push_back(mi.lmesh[vertxj*fullx+vertxi]);
    }

    for (int alpha=1; alpha<stencil_size; alpha++){
        int ydegree = alpha/polyn_order[0];
        int xdegree = alpha%polyn_order[0];

        for (int i = 0; i<stencil_size; i++){
            double * Dpq = CreateBasisPolynDeriv(xdegree,ydegree,i);

            for (int i_prime=0; i_prime<stencil_size; i_prime++){
                double * Dpq_prime = CreateBasisPolynDeriv(xdegree,ydegree,i_prime);

                for (int a=0; a<stencil_size; a++){ 
                for (int b=0; b<stencil_size; b++){
                    int ypow1 = a/polyn_order[0]; int xpow1 = a%polyn_order[0];
                    int ypow2 = b/polyn_order[0]; int xpow2 = b%polyn_order[0];

                    sigma[i_prime+i*stencil_size] += Dpq[a]*Dpq_prime[b]*
                                      NumIntegralFace(work,{xpow1+xpow2,ypow1+ypow2},stencil_center,h,poly);
                }}
                delete[] Dpq_prime;
            }
            delete[] Dpq;
        }
    }

}

double WenoStencil::ComputeSmoothnessIndicator(const MeshInfo& mi){

    double smoothIndi = 0.0;

    for (int pq=0; pq<stencil_size; pq++){
        valarray <int> cellpq = stencil_index_set[pq];
        cellpq[1] -= mi.ghost_vertx[1];
        cellpq[0] -= mi.ghost_vertx[0];

        for (int pq_prime=0; pq_prime<stencil_size; pq_prime++){
            valarray <int> cellpq_prime = stencil_index_set[pq_prime];
            cellpq_prime[1] -= mi.ghost_vertx[1];
            cellpq_prime[0] -= mi.ghost_vertx[0];

            smoothIndi += mi.localval[mi.localstart[1]+cellpq[1]][mi.localstart[0]+cellpq[0]]*
                          mi.localval[mi.localstart[1]+cellpq_prime[1]][mi.localstart[0]+cellpq_prime[0]]*
                          sigma[pq_prime+pq*stencil_size];
        }
    }

    return smoothIndi;
}

double WenoStencil::ComputeSmoothnessIndicator_Simple(const MeshInfo& mi){
    double smoothIndi = 0.0;

    for (auto & cell: stencil_index_set){
        valarray<int> shift = cell;
        shift[0] -= mi.ghost_vertx[0];
        shift[1] -= mi.ghost_vertx[1];
        smoothIndi += pow(mi.localval[mi.localstart[1]+targetCell_[1]][mi.localstart[1]+targetCell_[0]]- 
                                      mi.localval[mi.localstart[1]+shift[1]][mi.localstart[0]+shift[0]],2);
    }

    smoothIndi *= 1.0/(double)(stencil_size - 1);

    return smoothIndi;
}

/*
 *double WenoStencil::ComputeSmoothnessIndicator_Polyn(MeshInfo* mi){
 *
 *
 *
 *}
 */

void WenoStencil::CheckSigma(){
        for (int p=0; p<stencil_size; p++){
        for (int q=0; q<stencil_size; q++){
            cout << sigma[q+p*stencil_size] << "  ";
        }cout<<endl;}
}

void WenoStencil::PrintBasisPolyn(){
    int n = polyn_order[0]*polyn_order[1];
    for (int j=0; j<n; j++){
    for (int i=0; i<n; i++){
        printf("coeff = %.4f ",polyn[i+j*n]);
    }cout<<endl;}
}

// Weno Reconstruction ==============================================
WenoReconstruction::WenoReconstruction(const MeshInfo& mi, vector<double>& linWeights, vector<int *>& rangex, vector<int *>& rangey, point_index& target){

    assert(rangex.size()==rangey.size());
    assert(linWeights.size()==rangex.size());

    // Allocate required memories

    NonLinWeights_.resize(linWeights.size());
    fill(NonLinWeights_.begin(),NonLinWeights_.end(),0.0);

    omega_.resize(linWeights.size());
    fill(omega_.begin(),omega_.end(),0.0);

    sigma_.resize(linWeights.size());
    fill(sigma_.begin(),sigma_.end(),0.0);

    ws.resize(linWeights.size());

    double sumlinWeights=0.0;

    for (auto & w: linWeights){
        sumlinWeights += w;
    }

//    typedef WenoStencil * wsptr;

//    ws = new wsptr [linWeights.size()]; 

    for (int i=0; i<linWeights.size(); i++){
        linWeights_.push_back(linWeights[i]/sumlinWeights);
        // Create corresponding wenostencils
        ws[i]= new WenoStencil(mi,rangex[i],rangey[i],target);
    }

    rangex_ = rangex;
    rangey_ = rangey;

    target_ = target;
}

void WenoReconstruction::ComputeNonlinWeights(const MeshInfo& mi){

    for (int i=0; i<linWeights_.size(); i++){
        sigma_[i] = ws[i]->ComputeSmoothnessIndicator(mi);
        //sigma_[i] = ws[i]->ComputeSmoothnessIndicator_Simple(mi);
    }

    double sum_omega = 0.0;

    for (int i=0; i<linWeights_.size(); i++){
        int eta = max(ws[i]->polyn_order[0],ws[i]->polyn_order[1]);
        omega_[i] = linWeights_[i]/pow(sigma_[i]+ws[i]->Geth()/5.0,eta);
        sum_omega += omega_[i];
    }

    for (int i=0; i<linWeights_.size(); i++){
        NonLinWeights_[i] = omega_[i]/sum_omega;
    }

} 

double WenoReconstruction::WenoReconstStencil(const MeshInfo& mi, WenoStencil*& ws, point& target){

    double reconVal=0.0;

    for (int p=0; p<ws->stencil_size; p++){
        valarray<int> t = ws->stencil_index_set[p];
        t[0] = t[0] - mi.ghost_vertx[0];
        t[1] = t[1] - mi.ghost_vertx[1];
        for (int ypow=0; ypow<ws->polyn_order[1]; ypow++){
        for (int xpow=0; xpow<ws->polyn_order[0]; xpow++){
            int o = ypow*ws->polyn_order[0] + xpow;
            point shifted = (target-ws->GetCenter())/ws->Geth(); 
            //printf("(%d, %d), %3f ",t[0],t[1],mi.localval[t[1]][t[0]]);
            reconVal += mi.localval[mi.localstart[1]+t[1]][mi.localstart[0]+t[0]] * 
                        poly(shifted,{xpow,ypow}) * 
                        ws->polyn[o*ws->stencil_size+p]; 
        }}
    }

    return reconVal;
}

double WenoReconstruction::PointValueReconstruction(const MeshInfo& mi, point& target){

    double reconVal=0.0;

    for(int i=0; i<linWeights_.size(); i++){
        reconVal += NonLinWeights_[i]*WenoReconstStencil(mi, ws[i], target);
    }

    return reconVal; 
} 

double WenoReconstruction::Geth(){
    return ws[0]->Geth();
}

// Check calculated parameters
void WenoReconstruction::CheckSigma(){
    for (int i=0; i<linWeights_.size(); i++){

        ws[i]->CheckSigma();
        cout << endl;

    }cout << endl;
}

void WenoReconstruction::CheckSmoothnessIndicator(){
    cout << "Smooth Indicat:    ";
    for (int i=0; i<linWeights_.size(); i++){

        printf("%9f  ",sigma_[i]); 

    }cout << endl;
}

void WenoReconstruction::CheckPolynBasis(){
    for (int i=0; i<linWeights_.size(); i++){

        ws[i]->PrintBasisPolyn();
        cout << endl;

    }cout << endl;
}

void WenoReconstruction::CheckStencils(){
    for (int i=0; i<linWeights_.size(); i++){
        printf("[%2d,%2d]    " ,rangex_[i][0],rangex_[i][1]);
    } cout << endl;
    for (int i=0; i<linWeights_.size(); i++){
        printf("[%2d,%2d]    " ,rangey_[i][0],rangey_[i][1]);
    } cout << endl;

}

void WenoReconstruction::CheckNonlinWeights(){

    cout << "Nonlin Weights:    ";

    for (int i=0; i<linWeights_.size(); i++){

        printf("%9f  ",NonLinWeights_[i]); 

    } cout << endl; 
}

// ==================================================================

// 3D Case ==========================================================

