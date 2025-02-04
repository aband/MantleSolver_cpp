#include "advectiveflux.h"

inline double LFflux(const double& uin, const double& uout, 
                     const double& fin, const double& fout,
                     const double& LF){

    return 0.5*(fout + fin - LF*(uout - uin));
}

inline int derLFflux(const derivative& derin,  const derivative& derout,
                     const derivative& fderin, const derivative& fderout,
                     const double& LF,
                     derivative& work){

    // differentiating Lax_Friedrich flux

    work = fderin;
    unordered_map_arithmetic(work, fderout, std::plus<double>());

    unordered_map_arithmetic(work, derout, std::minus<double>(), 
                                   LF, std::multiplies<double>());

    unordered_map_arithmetic(work, derin,  std::plus<double>(), 
                                   LF, std::multiplies<double>());

    unordered_map_arithmetic(work, 0.5, std::multiplies<double>());

    return 1;
}

double edgefluxintegral(const MeshInfo& mi, 
                        const indice& gcellin,
                        const indice& gcellout,
                        const vertexSet& edge,
                        const Tensor<weights>& allwgts,
                        const vector<vertex>& vel,
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

    for (int g=0; g<gpe.size(); g++){
        vertex mapped = GaussMapPointsEdge({gpe[g]}, edge);
        double uin  = use.eval(mapped, ml, "all", 
                      allwgts({gcellin[0],gcellin[1]}), gcellin, lu); 
        double uout = use.eval(mapped, ml, "all", 
                      allwgts({gcellout[0],gcellout[1]}), gcellout, lu); 

        work += gwe[g] * LFflux(uin, uout, 
                                advfunc(uin,vel.at(g),unitNormal),
                                advfunc(uout,vel.at(g),unitNormal),1.0) * len/2.0; 
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

    const valarray<double>& gwe = GaussPointsEdge;
    vector<vertex> vel;
    vel.resize(gwe.size());
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
            flux    = edgefluxintegral(mi, globalcell, cellout, hori, allwgts, vel,
                                       ml, use, lu); 
        }

        horiedge({i,j}) = flux;
       
        vertexSet vert {corners.at(3), corners.at(0)};

        // boundary
        if (i==0){
            flux    = 0.0;
        } else {
            cellout = globalcell + mi.faceNormal[3];
            flux    = edgefluxintegral(mi, globalcell, cellout, vert, allwgts, vel,
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

// Similar way to obtain derivatives against u
int edgefluxintegral(const MeshInfo& mi, 
                     const indice& gcellin,
                     const indice& gcellout,
                     const vertexSet& edge,
                     const Tensor<weights>& allwgts,
                     const vector<vertex>& vel,
                     multilevel& ml,
                     mluse& use,
                     double ** lu,
                     derivative& der,
                     double& f){

    // Attention!!!!!! 
    // In serial code, stencil index equals global index
    // Which is not the case in parallel !!!!!!!!!!!
	 // fix it later

    f = 0.0;
    der.clear();

    //! Extract default gauess points and gauess weights.
    const valarray<double>& gwe = GaussWeightsEdge;
    const valarray<double>& gpe = GaussPointsEdge;

    // Get edge lendth and unit vector normal to the given edge
    double len = length(edge);
    vertex unitNormal = UnitNormal(edge,len);

    for (int g=0; g<gpe.size(); g++){

        vertex mapped = GaussMapPointsEdge({gpe[g]}, edge);

        derivative derin;
        use.der(mapped, ml, "all", allwgts({gcellin[0], gcellin[1]}), 
                gcellin, lu, mi, derin);

        derivative derout;
        use.der(mapped, ml, "all", allwgts({gcellout[0], gcellout[1]}), 
                gcellout, lu, mi, derout);

        double uin  = use.eval(mapped, ml, "all", 
                      allwgts({gcellin[0],gcellin[1]}), gcellin, lu); 
        double uout = use.eval(mapped, ml, "all", 
                      allwgts({gcellout[0],gcellout[1]}), gcellout, lu); 

        f += gwe[g] * LFflux(uin, uout, 
                             advfunc(uin,vel.at(g),unitNormal),
                             advfunc(uout,vel.at(g),unitNormal),1.0) * len/2.0; 
        derivative derfin;
        derivative derfout;
        derivative derLF;

        dadvfunc(derin , uin , vel.at(g), unitNormal, derfin);
        dadvfunc(derout, uout, vel.at(g), unitNormal, derfout);

//        cout << "derin : " << gcellin[0] << "  " << gcellin[1]<< endl;
//        unordered_map_print(derin);
//        cout << "derout : " << gcellout[0] << "  " << gcellout[1] << endl;
//        unordered_map_print(derout);
//        cout << "derfin : " << endl;
//        unordered_map_print(derfin);
//        cout << "derfout : " << endl;
//        unordered_map_print(derfout);

        derLFflux(derin, derout, derfin, derfout, 1.0, derLF);

        unordered_map_arithmetic(der, derLF, std::plus<double>(), 
                         gwe[g]*len/2.0, std::multiplies<double>());

    }

    return 1;
}

// Update edge flux and edge derivative of flux at once
int updateEdgeFlux(Tensor<double>& vertedge, Tensor<double>& horiedge,
                   Tensor<derivative>& vertedgeder, Tensor<derivative>& horiedgeder,
                   const MeshInfo& mi, double ** lu,
                   mluse& use, multilevel& ml, const Tensor<weights>& allwgts){

    const valarray<double>& gwe = GaussPointsEdge;
    vector<vertex> vel;
    vel.resize(gwe.size());
    for (int i=0; i<vel.size(); i++){
        vel.at(i) = {1,0};
    }

    // In this simplified test, no flux boundary are applied to all 
	 // four sides, flux = 0 and dflux = 0 as well
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
            edgefluxintegral(mi, globalcell, cellout, hori, 
            allwgts, vel, ml, use, lu, horiedgeder({i,j}), flux);
//            cout << i << "  " << j << endl;
//            unordered_map_print(horiedgeder({i,j})); cout << endl;
        }

        horiedge({i,j}) = flux;

        vertexSet vert {corners.at(3), corners.at(0)};

        // boundary
        if (i==0){
            flux    = 0.0;
        } else {
            cellout = globalcell + mi.faceNormal[3];
            edgefluxintegral(mi, globalcell, cellout, vert, allwgts, vel,
                             ml, use, lu, vertedgeder({i,j}), flux);

        }
        vertedge({i,j}) = flux;
    }}

    return 1;
}

int getcellflux(const MeshInfo& mi, const indice& gcell,
                const Tensor<double>& vertedge, 
                const Tensor<double>& horiedge,
                const Tensor<derivative>& vertedgeder,
                const Tensor<derivative>& horiedgeder,
                double& flux,
                derivative& dflux){

    double area = mi.cellArea.at(FlatIndic(mi,gcell));

    flux += horiedge({gcell[0], gcell[1]});

    flux -= horiedge({gcell[0], gcell[1]+1});

    flux += vertedge({gcell[0], gcell[1]});

    flux -= vertedge({gcell[0]+1, gcell[1]});

    flux /= area;

    // dflux 

    unordered_map_arithmetic(dflux, horiedgeder({gcell[0], gcell[1]}),
                             std::plus<double>());

    unordered_map_arithmetic(dflux, horiedgeder({gcell[0], gcell[1]+1}),
                             std::minus<double>());

    unordered_map_arithmetic(dflux, vertedgeder({gcell[0], gcell[1]}),
                             std::plus<double>());

    unordered_map_arithmetic(dflux, vertedgeder({gcell[0]+1, gcell[1]}),
                             std::minus<double>());

    unordered_map_arithmetic(dflux, 1.0/area,
                             std::multiplies<double>());

    return 1;
}
