#include "transport.h"

// Create Default (3,2) reconstruction
int TransportVariable::CreateDefaultReconstruction(const MeshInfo& mi){

	 int sizelgx = 3;
    int sizelgy = 3;
    int orderlg = 2;

	 int sizesmx = 2;
    int sizesmy = 2;
    int ordersm = 1;

    // (3,2) reconstruction but 1D
    vector<indice> sten_lg_pre = {{-1,-1}};
    vector<indice> sten_sm_pre = {{0,-1}, {0,0}, {-1,-1}, {-1,0}};

    CreateReconstruction(mi, sizelgx, sizelgy, orderlg, 
                             sizesmx, sizesmy, ordersm,
                             sten_lg_pre, sten_sm_pre); 

    return 1;
}

int TransportVariable::CreateReconstruction(const MeshInfo& mi, 
                                            int sizelgx, int sizelgy, int orderlg,
                                            int sizesmx, int sizesmy, int ordersm,
                                            vector<indice>& sten_lg_pre,
                                            vector<indice>& sten_sm_pre){

    int M = mi.MPIglobalCellSize[0];
    int N = mi.MPIglobalCellSize[1];

    int Mlg = M-sizelgx+1;
    int Nlg = N-sizelgy+1;

    int Msm = M-sizesmx+1;
    int Nsm = N-sizesmy+1;

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

    my_recon.resize(M*N);

    // Initializing reconstrucitons for each cell
	 // Change here to have sided reconstruction on the boundary
    for (int j=0; j<N; j++){
    for (int i=0; i<M; i++){
        int s = j*M+i;
 
        my_recon.at(s) = new reconstruction();

        my_recon.at(s)->init(sizesmx,sizesmy,
                             sizelgx,sizelgy,
                             ordersm,orderlg,
                             sten_lg_pre, sten_sm_pre, mi,{i,j});
    }}

    return 1;
}

// Update reconstruction nonlinear weights and smoothness indicators
int TransportVariable::UpdateRecon(const MeshInfo& mi, double ** locvals){

    int M = mi.MPIglobalCellSize[0];
    int N = mi.MPIglobalCellSize[1];

    sigma_lg.clear();
    sigma_lg.resize(stenlg.size());

    sigma_sm.clear();
    sigma_sm.resize(stensm.size());

    for (int s=0; s<stenlg.size(); s++){
        sigma_lg.at(s) = stenlg.at(s).sigma(locvals);
    }

    for (int s=0; s<stensm.size(); s++){
       sigma_sm.at(s) = stensm.at(s).sigma(locvals);
    }

    // Setup nonlinear weights
    for (int s=0; s<my_recon.size(); s++){
        my_recon.at(s)->extractsigma(sigma_lg, sigma_sm);
        my_recon.at(s)->setWgts(1.0/(double)M/(double)N);
    }

    return 1;
}

// Evaluate at each gauss points on the edges
int TransportVariable::Evaluate(const MeshInfo& mi, DM dmu){

    Vec localu;
    DMGetLocalVector(dmu, &localu);

    DMGlobalToLocalBegin(dmu, sol, INSERT_VALUES, localu);
    DMGlobalToLocalEnd(dmu, sol, INSERT_VALUES, localu); 

    double ** lu;
    DMDAVecGetArray(dmu, localu, &lu);

    UpdateRecon(mi, lu);

    for (int j=0; j<mi.MPIglobalCellSize[1]; j++){
    for (int i=0; i<mi.MPIglobalCellSize[0]; i++){

        EvaluateEdge(i, j, mi, lu);

    }}

    DMDAVecRestoreArray(dmu, localu, &lu);
    DMRestoreLocalVector(dmu, &localu);

    return 1;
}

// Evaluate reconstruction value at each gauss points on the edges
// With the same order left - bottom - right - top
int TransportVariable::EvaluateEdge(int i, int j, const MeshInfo& mi, double ** locvals){

    // Get gauss points on the edges 
    const valarray<double>& gwe = GaussWeightsEdge;
    const valarray<double>& gpe = GaussPointsEdge;

    vertexSet corners = extractCorners(mi, {i,j});

    vector<vertex> gaussp;
    gaussp.resize(4*gpe.size());

    vector<indice> edgeBound {{3,0},{0,1},{2,1},{3,2}};

    for (int e=0; e<4; e++){
        vertexSet edge = {corners.at(edgeBound.at(e)[0]), corners.at(edgeBound.at(e)[1])};

        for (int g=0; g<gpe.size(); g++){
            gaussp.at(e*gpe.size()+g) = GaussMapPointsEdge({gpe[g]}, edge);
        }
    }

    my_recon.at(j*mi.MPIglobalCellSize[0]+i)->eval(locvals, gaussp, stenlg, stensm); 

    return 1;
}
