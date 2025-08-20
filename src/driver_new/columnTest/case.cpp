#include "driver.h"
// Initialize dimensionless composition and enthalpy
double InitCD(const valarray<double>& point,
              const vector<double>& param){

    // Constant composition value 
    return 0.04;
}

double InitHD(const valarray<double>& point,
              const vector<double>& param){

    // Linear simple distribution of enthalpy
	 // We pass nondimensionalize normalization factor in param.at(0)
    //double HD = 0.01;

    double HD = 2.9-2.5*point[1];

    if (point[1] < -0.20){HD = 2.9 + 2.5*0.20;}

    return HD;
}

// Case specified different case requires different 
// Prepare reconstruction stencils
int Driver::PrepareTransport(double (*initHD)(const valarray<double>& point, 
                                              const vector<double>& param),
                             double (*initCD)(const valarray<double>& point, 
                                              const vector<double>& param)){

    PetscCall(DMCreateGlobalVector(dmu, &globalCD));
    PetscCall(DMCreateGlobalVector(dmu, &globalHD));

    // Assign cell averaged values as initial condition
    SimpleInitialValue(dmMesh, dmu, &globalmesh, &globalCD, {H_,0.0}, initCD);
    SimpleInitialValue(dmMesh, dmu, &globalmesh, &globalHD, {myPhase->pp->l0*H_,-0.7*H_}, initHD);

    int M = mi.MPIglobalCellSize[0];
    int N = mi.MPIglobalCellSize[1];

    stenlg.resize((M-2)*(N-2));

    for (int j=0; j<N-2; j++){
    for (int i=0; i<M-2; i++){
        int s = j*(M-2)+i;
        stenlg.at(s) = tensorstencilpoly(2, 3, 3);
        stenlg.at(s).setCoef(mi,i,j);
		  stenlg.at(s).setSigma();
		  stenlg.at(s).startx = i;
		  stenlg.at(s).starty = j;
    }}
  
    stensm.resize((M-1)*(N-1));

    for (int j=0; j<N-1; j++){
    for (int i=0; i<M-1; i++){
        int s = j*(M-1) + i;
        stensm.at(s) = tensorstencilpoly(1, 2, 2);
        stensm.at(s).setCoef(mi,i,j);
		  stensm.at(s).setSigma();
		  stensm.at(s).startx = i;
		  stensm.at(s).starty = j;
    }}

    // (3,2) reconstruction but 1D
    //vector<indice> sten_lg_pre = {{0,-1}};
    //vector<indice> sten_sm_pre = {{0,-1}, {0,0}};

    vector<indice> sten_lg_pre = {{-1,-1}};
    vector<indice> sten_sm_pre = {{-1,-1}, {0,0}, {-1,0}, {0,-1}};

    my_recon_HD.resize(M*N);
    my_recon_CD.resize(M*N);

    // Initializing reconstrucitons for each cell
    for (int j=0; j<N; j++){
    for (int i=0; i<M; i++){
        int s = j*M+i;
 
        my_recon_HD.at(s) = reconstruction();

        my_recon_HD.at(s).init(2,2,3,3,1,2,sten_lg_pre, sten_sm_pre, mi,{i,j});

        my_recon_CD.at(s) = reconstruction();

        my_recon_CD.at(s).init(2,2,3,3,1,2,sten_lg_pre, sten_sm_pre, mi,{i,j});
    }}

    // Compute bottom fixed value
    HDbottom = initHD({0.0,-1*H_},{myPhase->pp->l0*H_,-0.7*H_});
    CDbottom = 0.04;

    return 1;
}
