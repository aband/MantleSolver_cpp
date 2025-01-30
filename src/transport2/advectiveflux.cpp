#include "advectiveflux.h"

inline double LFflux(const double& uin,  const double& uout, 
                     const double& fin, const double& fout,
                     const double& LF){

    return 0.5*(fout + fin - LF*(uout - uin));
}

double edgefluxintegral(const MeshInfo& mi, 
                        const indice& gcellin,
                        const indice& gcellout,
                        const vertexSet& edge,
                        const Tensor<weights>& allwgts,
                        multilevel& ml,
                        mluse& use,
                        double ** lu){

    double work = 0.0;

    //! Extract default gauess points and gauess weights.
    const valarray<double>& gwe = GaussWeightsEdge;
    const valarray<double>& gpe = GaussPointsEdge;

    // Get edge lendth and unit vector normal to the given edge
    double len = length(edge);
    vertex unitNormal = UnitNormal(edge,len);

    // Constant velocity for testing
    vertex vel = {1,0};

    for (int g=0; g<gpe.size(); g++){
        vertex mapped = GaussMapPointsEdge({gpe[g]}, edge);
        double uin  = use.eval(mapped, ml, "all", 
                      allwgts({gcellin[0],gcellin[1]}), gcellin, lu); 
        double uout = use.eval(mapped, ml, "all", 
                      allwgts({gcellout[0],gcellout[1]}), gcellout, lu); 

        work += gwe[g] * LFflux(uin, uout, 
                                advfunc(uin,vel,unitNormal),
                                advfunc(uout,vel,unitNormal),1.0) * len/2.0; 
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
            flux    = edgefluxintegral(mi, globalcell, cellout, hori, allwgts,
                                       ml, use, lu); 
        }

        horiedge({i,j}) = flux;
       
        vertexSet vert {corners.at(3), corners.at(0)};

        // boundary
        if (i==0){
            flux    = 0.0;
        } else {
            cellout = globalcell + mi.faceNormal[3];
            flux    = edgefluxintegral(mi, globalcell, cellout, vert, allwgts, 
                                       ml, use, lu);
        }
        vertedge({i,j}) = flux;
    }}

    return 1;
}

double getcellflux(const MeshInfo& mi, const indice& gcell,
                   const Tensor<double>& vertedge, 
                   const Tensor<double>& horiedge){

    double work = 0.0;

    double area = mi.cellArea.at(FlatIndic(mi,gcell));

    work += horiedge({gcell[0], gcell[1]});

    work -= horiedge({gcell[0], gcell[1]+1});

    work += vertedge({gcell[0], gcell[1]});

    work -= vertedge({gcell[0]+1, gcell[1]});

    work /= area;

    //cout << "At cell " << gcell[0] << "  " << gcell[1] << endl;
    //cout << "left   : " << horiedge({gcell[0], gcell[1]}) << " ";
    //cout << "right  : " << horiedge({gcell[0], gcell[1]+1}) << " ";
    //cout << "bottom : " << vertedge({gcell[0], gcell[1]}) << " ";
    //cout << "top    : " << vertedge({gcell[0]+1, gcell[1]}) << " ";
    //cout << endl << endl;;

    return work;  
}
