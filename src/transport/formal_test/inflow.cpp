#include "advectiveflux.h"
#include "trans_param.h"
#include "temp.h"

inline double tempLFflux(const double& uin, const double& uout, 
                         const double& fin, const double& fout,
                         const double& LF){

    return 0.5*(fout + fin - LF*(uout - uin));
}

double edgefluxintegral_inflow(const MeshInfo& mi, 
                               const indice& gcellin,
                               const indice& gcellout,
                               const vertexSet& edge,
                               const Tensor<weights>& allwgts,
                               const vector<vertex>& vel,
                               multilevel& ml,
                               mluse& use,
                               double ** lu){

    double work = 0.0;

    //! Extract 
    const valarray<double>& gwe = GaussWeightsEdge;
    const valarray<double>& gpe = GaussPointsEdge;

	 double len = length(edge);

    vertex unitNormal = UnitNormal(edge,len);

    for (int g=0; g<gpe.size(); g++){
        vertex mapped = GaussMapPointsEdge({gpe[g]}, edge);
        double uin  = use.eval(mapped, ml, location(mi,gcellin), 
                      allwgts({gcellin[0],gcellin[1]}), gcellin, lu); 
        double uout = use.eval(mapped, ml, location(mi,gcellout), 
                      allwgts({gcellout[0],gcellout[1]}), gcellout, lu); 
        //double LF = sqrt(vel.at(g)[0]*vel.at(g)[0] + vel.at(g)[1]*vel.at(g)[1]);
        double LF = abs(vel.at(g)[0]*unitNormal[0] + vel.at(g)[1]*unitNormal[1]);
        //LF = find_max(abs(dfdu(uin)), abs(dfdu(uout))) * LF;
        LF = 1.0;

        work += gwe[g] * tempLFflux(uin, uout, 
                                    advfunc(uin,vel.at(g),unitNormal),
                                    advfunc(uout,vel.at(g),unitNormal),LF) 
                       * len/2.0; 
    }

    return work;
}

double func(const vertex& point,
            const vector<double>& param){

//    if (point[1] < 0.25 || point[1] > 0.75){
//        return pow(sin(2*M_PI*point[1]),2);
//    } else {
//        return 1.0;
//    }

    return 0.0;
}

double dfdu(const double& u){

    return u*u*u;
}

// Burgers for testing
double advfunc(const double& u, 
               const vertex& vel, const vertex& unitnormal){

    // A Burgers type flux
    return u*u*u*u/4.0 *(unitnormal[0]*vel[0] + unitnormal[1]*vel[1]);
}

int dadvfunc(const derivative& du, const double& u, const vertex& vel, const vertex& unitnormal, derivative& work){

    // compute df/du = df/dR * dR/du

    double direction = dfdu(u)*(unitnormal[0]*vel[0] + unitnormal[1]*vel[1]);

    unordered_map_arithmetic(work, du, std::plus<double>(), direction, std::multiplies<double>());

    return 1;
}

double inflowfunc(const vertex& point, const vector<double>& param){

    double work=0.0;

    if (point[1] < 0.25 || point[1] > 0.75){
        return -pow(sin(2*M_PI*point[1]),2);
    } else {
        return -1.0;
    }

    return work;
}

int updateEdgeFlux(Tensor<double>& vertedge, Tensor<double>& horiedge,
                   const MeshInfo& mi, double ** lu,
                   mluse& use, multilevel& ml, const Tensor<weights>& allwgts){

    // Udpate every left and bottom edge for each cell
    //for (int j=mi.MPIlocalCellStart[1]; j<mi.MPIlocalCellStart[1] + mi.MPIlocalCellSize[1] + 1; j++){
    //for (int i=mi.MPIlocalCellStart[0]; i<mi.MPIlocalCellStart[0] + mi.MPIlocalCellSize[0] + 1; i++){
	 // Test for serial now
    Tensor_zero(vertedge);
    Tensor_zero(horiedge);
   
    const valarray<double>& gpe = GaussPointsEdge;
    vector<vertex> vel;
    vel.resize(gpe.size());
    for (int i=0; i<vel.size(); i++){
        vel.at(i) = {1,0};
    }

    for (int j=0; j<mi.MPIglobalCellSize[1]; j++){
    for (int i=0; i<mi.MPIglobalCellSize[0]; i++){

        double flux = 0.0;
        indice globalcell {i,j};
        indice cellout;

        // Extract corners with respect to given global indice
        vertexSet corners = extractCorners(mi, globalcell); 

        // Compute and restore Horizontal flux
        vertexSet hori {corners.at(0), corners.at(1)};
        // boundary
        if (j==0){
            // Temperatory
            flux    = 0.0;
        } else {
            cellout = globalcell + mi.faceNormal[0];
            //flux    = edgefluxintegral_inflow(mi, globalcell, cellout,hori,allwgts,vel,
            //                           ml, use, lu);
				flux = 0.0;
        }

        horiedge({i,j}) = flux;
       
        vertexSet vert {corners.at(3), corners.at(0)};

        // boundary
        if (i==0){
            flux = 0.0;

        } else {
            cellout = globalcell + mi.faceNormal[3];
            flux    = edgefluxintegral_inflow(mi, globalcell, cellout, vert, allwgts, vel,
                                            ml, use, lu);
        }
        vertedge({i,j}) = flux;
    }}

    // Free outflow
    for (int j=0; j<mi.MPIglobalCellSize[1]; j++){

        indice gcell {mi.MPIglobalCellSize[0]-1, j};
        vertexSet corners = extractCorners(mi, gcell);
        vertexSet vert    = {corners.at(2), corners.at(1)};

        vertedge({mi.MPIglobalCellSize[0], j}) = edgefluxintegral(mi, gcell, vert, allwgts, vel, ml, use, lu);
    }

    // Prescribed inflow boundary
    for (int j=0; j<mi.MPIglobalCellSize[1]; j++){

        indice gcell {0, j};
        vertexSet corners = extractCorners(mi, gcell);
        vertexSet vert    = {corners.at(3), corners.at(0)};

        vertedge({0,j}) = edgefluxintegral(vert, inflowfunc, {0}, vel);
    }

    return 1;
}
