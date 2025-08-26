#include "driver.h"
double dfdu(const double& u){

    return 1.0;
}

// Burgers for testing
double advfunc(const double& u, 
               const vertex& vel, const vertex& unitnormal){

    // A Burgers type flux
//    return u*u/2.0 *(unitnormal[0]*vel[0] + unitnormal[1]*vel[1]);
    return u *(unitnormal[0]*vel[0] + unitnormal[1]*vel[1]);
}

int dadvfunc(const derivative& du, const double& u, const vertex& vel, const vertex& unitnormal, derivative& work){

    // compute df/du = df/dR * dR/du

    double direction = dfdu(u)*(unitnormal[0]*vel[0] + unitnormal[1]*vel[1]);

    unordered_map_arithmetic(work, du, std::plus<double>(), direction, std::multiplies<double>());

    return 1;
}

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

double inflow(const valarray<double>& point,
              const vector<double>& param){

    return param[0];
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

	 int sizelgx = 1;
    int sizelgy = 3;
    int orderlg = 2;
    int Mlg = M-sizelgx+1;
    int Nlg = N-sizelgy+1;

	 int sizesmx = 1;
    int sizesmy = 2;
    int ordersm = 1;
    int Msm = M-sizesmx+1;
    int Nsm = N-sizesmy+1;

    // (3,2) reconstruction but 1D
    vector<indice> sten_lg_pre = {{0,-1}};
    vector<indice> sten_sm_pre = {{0,-1}, {0,0}};

    //vector<indice> sten_lg_pre = {{-1,-1}};
    //vector<indice> sten_sm_pre = {{-1,-1}, {0,0}, {-1,0}, {0,-1}};

    // ==================================================
    stenlg.resize(Mlg*Nlg);

    for (int j=0; j<Nlg; j++){
    for (int i=0; i<Mlg; i++){
        int s = j*Mlg+i;
        stenlg.at(s) = tensorstencilpoly(orderlg, sizelgx, sizelgy);
        stenlg.at(s).setCoef(mi,i,j);
		  stenlg.at(s).setSigma();
		  stenlg.at(s).startx = i;
		  stenlg.at(s).starty = j;
    }}
  
    stensm.resize(Msm*Nsm);

    for (int j=0; j<Nsm; j++){
    for (int i=0; i<Msm; i++){
        int s = j*Msm + i;
        stensm.at(s) = tensorstencilpoly(ordersm, sizesmx, sizesmy);
        stensm.at(s).setCoef(mi,i,j);
		  stensm.at(s).setSigma();
		  stensm.at(s).startx = i;
		  stensm.at(s).starty = j;
    }}

    my_recon_HD.resize(M*N);
    my_recon_CD.resize(M*N);

    // Initializing reconstrucitons for each cell
    for (int j=0; j<N; j++){
    for (int i=0; i<M; i++){
        int s = j*M+i;
 
        my_recon_HD.at(s) = reconstruction();
        my_recon_CD.at(s) = reconstruction();

        my_recon_HD.at(s).init(sizesmx,sizesmy,sizelgx,sizelgy,ordersm,orderlg,sten_lg_pre, sten_sm_pre, mi,{i,j});
        my_recon_CD.at(s).init(sizesmx,sizesmy,sizelgx,sizelgy,ordersm,orderlg,sten_lg_pre, sten_sm_pre, mi,{i,j});
    }}

    // Compute bottom fixed value
    HDbottom = initHD({0.0,-1*H_},{myPhase->pp->l0*H_,-0.7*H_});
    CDbottom = 0.04;

    return 1;
}

int Driver::computeEdgeFlux(vector<double>& edgefluxHD, 
                            vector<double>& edgefluxCD, 
                            double t, 
                            double ** lHD, double ** lCD){

    int M = mi.MPIglobalCellSize[0];
    int N = mi.MPIglobalCellSize[1];

    edgefluxHD.clear();
    edgefluxHD.resize(M*(N+1) + N*(M+1));
    std::fill(edgefluxHD.begin(), edgefluxHD.end(), 0);

    edgefluxCD.clear();
    edgefluxCD.resize(M*(N+1) + N*(M+1));
    std::fill(edgefluxCD.begin(), edgefluxCD.end(), 0);

    const valarray<double>& gwe = GaussWeightsEdge;
    const valarray<double>& gpe = GaussPointsEdge;

    std::vector<vertex> gaussp;
    gaussp.resize(gpe.size());

    // Effective velocity should be used for transport of concentration
    vector<vertex> effvel;   effvel.resize(gaussp.size());
    vector<vertex> phasevel; phasevel.resize(gaussp.size());
    vector<vertex> solidvel; phasevel.resize(gaussp.size());
    vector<double> TDin;    TDin.resize(gaussp.size());
    vector<double> TDout;   TDout.resize(gaussp.size());
    vector<double> dTdHin;  dTdHin.resize(gaussp.size());
    vector<double> dTdHout; dTdHout.resize(gaussp.size());

    vector<double> CDin;  CDin.resize(gaussp.size());
    vector<double> CDout; CDout.resize(gaussp.size());
    vector<double> HDin;  HDin.resize(gaussp.size());
    vector<double> HDout; HDout.resize(gaussp.size());

    vector<double> nonlinuin;     nonlinuin.resize(gaussp.size());
    vector<double> nonlinuout;    nonlinuout.resize(gaussp.size());
    vector<double> nonlinfin;     nonlinfin.resize(gaussp.size());
    vector<double> nonlinfout;    nonlinfout.resize(gaussp.size());
    vector<double> nonlindfduin;  nonlindfduin.resize(gaussp.size());
    vector<double> nonlindfduout; nonlindfduout.resize(gaussp.size());
    vector<vertex> nonlinvel;     nonlinvel.resize(gaussp.size());

    // 1D column test, asserting flux across the vertical edge to be zero
    for (int j=0; j<N-1; j++) {
    for (int i=0; i<M  ; i++) {

        clearall(effvel, phasevel, solidvel, TDin, TDout, dTdHin, dTdHout, 
                 CDin, CDout, HDin, HDout);

        vertexSet corners = extractCorners(mi, {i,j});
        vertexSet edge {corners.at(3), corners.at(2)};

        // extract velocity on this edge 
        for (int g=0; g<gpe.size(); g++){
            gaussp.at(g) = GaussMapPointsEdge({gpe[g]},edge);
        }   

        // Negative cell and positive cell
        int neg = j*M + i; 
        int pos = (j+1)*M + i;

        int position = (j+1)*M +i;

        computeEffVel(gaussp, edge, {i,j}, {i,j+1}, lHD, 
                      lCD, effvel, phasevel, solidvel, 
                      TDin, TDout, dTdHin, dTdHout, CDin, CDout, HDin, HDout);

        edgefluxHD.at(position) = advflux_edge(my_recon_HD.at(neg), my_recon_HD.at(pos), stenlg, stensm, lHD, edge, phasevel, true, 1.0);

        edgefluxCD.at(position) = advflux_edge(my_recon_CD.at(neg), my_recon_CD.at(pos), stenlg, stensm, lCD, edge, effvel, true, 1.0);

    }}

    // Bottom is dirichlet boundary condition
    // Top is free outflow boundary condition
    for (int i=0; i<M; i++){
        vertexSet corners = extractCorners(mi, {i, 0});
        vertexSet edge {corners.at(0), corners.at(1)};

        int position = i;

        clearall(effvel, phasevel, solidvel, TDin, TDout, dTdHin, dTdHout, 
                 CDin, CDout, HDin, HDout);

        // extract velocity on this edge 
        for (int g=0; g<gpe.size(); g++){
            gaussp.at(g) = GaussMapPointsEdge({gpe[g]},edge);
        }   

        computeEffVel(gaussp, edge, {i,0}, lHD, lCD, 
                      effvel, phasevel, solidvel, TDin, dTdHin, 
                      CDin, HDin);

        edgefluxHD.at(position) = advflux_edge(inflow, {HDbottom}, edge, phasevel, true, 1.0);
        edgefluxCD.at(position) = advflux_edge(inflow, {CDbottom}, edge, effvel,   true, 1.0);

    }

    for (int i=0; i<M; i++){

        vertexSet corners = extractCorners(mi, {i, N-1});
        vertexSet edge {corners.at(3), corners.at(1)};

        int position = N*M + i;
 
        int cell = (N-1)*M + i;

        clearall(effvel, phasevel, solidvel, TDin, TDout, dTdHin, dTdHout, 
                 CDin, CDout, HDin, HDout);

        // extract velocity on this edge 
        for (int g=0; g<gpe.size(); g++){
            gaussp.at(g) = GaussMapPointsEdge({gpe[g]},edge);
        }   
   
        computeEffVel(gaussp, edge, {i,N-1}, lHD, lCD, 
                      effvel, phasevel, solidvel, TDin, dTdHin, 
                      CDin, HDin);

        edgefluxHD.at(position) = advflux_edge(my_recon_HD.at(cell), stenlg, stensm, lHD, edge, phasevel, true, 1.0);
        edgefluxCD.at(position) = advflux_edge(my_recon_CD.at(cell), stenlg, stensm, lCD, edge, effvel,   true, 1.0);
    }

    return 1;
}
