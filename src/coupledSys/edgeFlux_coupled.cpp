#include "edgeFlux.h"
#include "preconst.h"

inline bool OutBndryCell(const MeshInfo& mi, 
                         const indice& gcell){

    // Check if the Cell is out of domain or not.
    bool work = false;

    if (gcell[0] < 0 || gcell[0] > mi.MPIglobalCellSize[0]-1 || 
        gcell[1] < 0 || gcell[1] > mi.MPIglobalCellSize[1]-1){ 

        work = true;
    }

    return work;
}

// pick the cell index that inside the compuitational domain
inline indice PickCellInside(const MeshInfo& mi,
                             const indice& gCellIn,
                             const indice& gCellOut){

    if (OutBndryCell(mi, gCellIn)){
        return gCellOut;
    } else if (OutBndryCell(mi, gCellOut)){
        return gCellIn;
    } else {
        // Both gCellIn and gCellOut are inside boundary
        return gCellIn;
    }

}

int extractEdgeGaussPoints(unordered_map<int, vector<double>>& edgeGaussPointsAll,
                           const MeshInfo& mi,
                           const std::valarray<double>& gpe){

}

int extractVelocityAll(unordered_map<int, vector<vertex>>& velocityAll,
                       const unordered_map<int, vector<double>>& edgeGaussPointsAll,
                       Vec * vel, Vec * g){


}

inline double computeFlux(const MeshInfo& mi,
                          const MLWENO::MLWENOUse& mlu,
                          const indice& globalCellIn,
                          const indice& globalCellOut,
                          const std::string& locationIn,
                          const std::string& locationOut,
                          const std::vector<vertex>& edge,
                          const std::vector<vertex>& gauss_p,
                          const std::vector<vertex>& vel_darcy,
                          const std::vector<vertex>& vel_stokes,
                          const std::valarray<double>& gwe,
                          const std::valarray<double>& gpe,
                          const int& edgemark,
                          fluxFunc       fluxfunc,
                          fluxFuncBndry  fluxfuncbndry){
    // Compute integral of flux on edge
    // mluIn used for value where unit normal points to
    // mluOut used for value where unit normal tails at

    // Calculate unit normal vector and length of corresponding edge
    double len = length(edge);
    vertex unitNormal = UnitNormal(edge, len); 

    double flux = 0.0;
    vector<double> fluxpoint;
    vector<double> directionIn;
    vector<double> directionOut;
    vector<double> LFParam;

    int flag = 0;

    for (const auto& it : velocityIn){
        directionIn.push_back(it[0]*unitNormal[0] + 
                              it[1]*unitNormal[1]); 
    }

    for (const auto& it : velocityOut){
        directionOut.push_back(it[0]*unitNormal[0] + 
                               it[1]*unitNormal[1]);
    }

    // Implement local LF scheme here
    // max |df/du| will be used as LF stabilizer
    LFParam.resize(fluxpoint.size());

    for (int i=0; i<fluxpoint.size(); i++){
        LFParam.at(i) = find_max<double>(directionIn.at(i), directionOut.at(i));
    }

    // Distinguish different boundary conditions
    if(OutBndryCell(mi, globalCellIn)){
        flag = 1;
        fluxpoint = fluxfuncbndry(mi, mlu, edge, unitNormal, len, globalCellOut,
                                 locationOut, directionOut, directionOut, gpe, 
                                 AssignBoundary(globalCellOut, locedge), flag, locedge) ;

    }else if (OutBndryCell(mi, globalCellOut)){
        flag = 0;
        fluxpoint = fluxfuncbndry(mi, mlu, edge, unitNormal, len, globalCellIn,
                                  locationIn, directionIn, directionIn, gpe, 
                                  AssignBoundary(globalCellIn, locedge), flag, locedge) ;

    } else {
        // Interior edge
        fluxpoint = fluxfunc(mi, mlu, edge, unitNormal, 
                             len, globalCellIn, globalCellOut,
                             locationIn, locationOut, LFParam, directionIn, directionOut, gpe) ;
    }

    // Integrate with gauess quadratrue scheme
    for (int g=0; g<gwe.size(); g++){
        flux += gwe[g] * fluxpoint.at(g); 
    }

    return flux;
}

inline int extractVertEdgeInfo(const MeshInfo& mi, 
                               const indice& local,
                               const indice& ghostShift,
                               indice& gCellOut,
                               indice& gCellIn,
                               edgeEnds<vertex> edgeEndsVertex,
                               edgeEnds<indice> edgeEndsIndice){

    // Extract information for vertical edges
    // index counted from top to bottom
    // Cell on right of the edge is regarded as "In" Cell
    // Cell on left of the edge is regarded as "Out" Cell
    gCellIn  = local + mi.MPIlocalCellStart; 
    gCellOut = {gCellIn[0] - 1, gCellIn[0]};

    edgeEndsIndice.end   = local + ghostShift;
    edgeEndsIndice.start = {edgeEndsIndice.end[0], edgeEndsIndice.end[1]-1};

    edgeEndsVertex.start = mi.lmesh[FlatIndic(mi.MPIlocalVertexSizeFull[0], edgeEndsIndice.start)];
    edgeEndsVertex.end   = mi.lmesh[FlatIndic(mi.MPIlocalVertexSizeFull[0], edgeEndsIndice.end)];

    // 1 represents vertical edge
    return 1;
}

inline int extractHoriEdgeInfo(const MeshInfo& mi,
                               const indice& local,
                               const indice& ghostShift,
                               indice& gCellOut,
                               indice& gCellIn,
                               edgeEnds<vertex>& edgeEndsVertex,
                               edgeEnds<indice>& edgeEndsIndice){

    // Extract information for horizontal edges
    // index counted from left to right
    // Cell on top of the edge is regarded as "In" cell 
    // Cell on bottom of the edge is regarded as "Out" cell
    // Unit normal vector pointing from top to bottom

    gCellIn = local + mi.MPIlocalCellStart;
    gCellOut= {gCellIn[0],gCellIn[1] - 1};

    edgeEndsIndice.start = local + ghostShift; 
    edgeEndsIndice.end   = {edgeEndsIndice.start[0] + 1, edgeEndsIndice.start[0]};

    edgeEndsVertex.start = mi.lmesh[FlatIndic(mi.MPIlocalVertexSizeFull[0], edgeEndsIndice.start)];
    edgeEndsVertex.end   = mi.lmesh[FlatIndic(mi.MPIlocalVertexSizeFull[0], edgeEndsIndice.end)];

    // 2 represents horizontal edge
    return 2;
}

typedef int (*extractEdgeInfoFunc) (const MeshInfo& mi,
                                    const indice& local,
                                    const indice& ghostShift,
                                    indice& gCellOut,
                                    indice& gCellIn,
                                    edgeEnds<vertex>& edgeEndsVertex,
                                    edgeEnds<indice>& edgeEndsIndice);

inline int getEdgeFlux(const MeshInfo& mi,
                       const int& i,
                       const int& j,
                       double * edgeFlux,
                       const int& shift,
                       const MLWENO::MLWENOUse& mlu,
                       indice& gCellOut,
                       indice& gCellIn,
                       const valarray<double>& gwe,
                       const valarray<double>& gpe,
                       edgeEnds<vertex>& edgeEndsVertex,
                       edgeEnds<indice>& edgeEndsIndice,
                       vector<vertex>& gauss_p,
                       fluxFunc      fluxfunc,
                       fluxFuncBndry fluxfuncbndry,
                       extractEdgeInfoFunc extractEdgeInfo,
                       Phase * phase,
                       Vec * sol_darcy, Vec * g_darcy,
                       Vec * sol_stokes, Vec * g_stokes,
                       const int* refmap_darcy,
                       const int* refmap_stokes,
                       basis& mybasis,
                       BRMixed& br,
                       Hdivmixed& hdiv){

    indice ghostShift {mi.vertexGhostLayerSize, mi.vertexGhostLayerSize};

    // Compute flux on the edge
    indice local {i,j};
    int flatlocal = FlatIndic(mi.MPIlocalCellSize[0], local);

    // Marks whether it is a horizontal or vertical edge
    int edgemark = extractEdgeInfo(mi, local, ghostShift, gCellOut, gCellIn, edgeEndsVertex, edgeEndsIndice);

    vector<vertex> edge {edgeEndsVertex.start, edgeEndsVertex.end};

    double len = length(edge);
    vertex unitNormal = UnitNormal(edge, len);

    // Get gauss quadrature points
    for (int g=0; g<gpe.size(); g++){gauss_p.at(g) = GaussMapPointsEdge({gpe[g]}, edge);}

    // Pick inside cell and flag
    indice gcell_inside = PickCellInside(mi, gCellIn, gCellOut);

    // Extract velocity
    // unscaled darcy
    vector<vertex> vel_darcy  = ExtractVelocity(sol_darcy, g_darcy, refmap_darcy, 
                                      mi, gauss_p, gcell_inside, hdiv, mybasis);

    vector<vertex> vel_stokes = ExtractVelocity(sol_stokes, g_stokes, refmap_stokes, 
                                      mi, gauss_p, gcell_inside, br, mybasis);

    // Compute edge flux

    edgeFlux[flatlocal + shift] = computeFlux(mi, mlu, globalCellIn, gobalCellOut, 
                                              locationIn, locationOut, edge, gauss_p, 
                                              vel_darcy, vel_stokes, gwe, gpe, edgemark, fluxfunc, fluxfuncbndry);

    return 1;
}

double * edgeFluxAll(const MeshInfo& mi,
                     Vec * sol_darcy, Vec * g_darcy,
                     Vec * sol_stokes, Vec * g_stokes,
                     const int* refmap_darcy,
                     const int* refmap_stokes,
                     const MLWENO::MLWENOUse& mlu,
                     basis& mybasis,
                     BRMixed& br,
                     Hdivmixed& hdiv,
                     Phase * phase,
                     fluxFunc      fluxfunc,
                     fluxFuncBndry fluxfuncbndry){

    // globalcell in and global cell out are fixed with the point direction of
    // the normal vector on the edge.
    // For vertical edges, left  cell is labeled in,
    //                     right cell is labeled out
	 //                     e.g. unit normal pointing from left to right
    // For horizontal edges, up cell is labeled in,
	 //                       bottom cell is labeled out
	 //                     e.g. unit normal pointing from up to bottom

    double * edgeflux = (double *)malloc(sizeof(double) * 
                        mi.MPIlocalHoriEdgeSize * mi.MPIlocalVertEdgeSize);

    vertex start, end;

    indice ghostShift {mi.vertexGhostLayerSize, mi.vertexGhostLayerSize};

    // Get gauss points and weights
    const valarray<double>& gwe = GaussWeightsEdge;
    const valarray<double>& gpe = GaussPointsEdge;

    vector<vertex> gauss_p; 
    gauss_p.resize(gpe.size());

    indice gCellOut;
    indice gCellIn;
    edgeEnds<vertex> edgeEndsVertex;
    edgeEnds<indice> edgeEndsIndice;

    // Loop through the entire local mesh chunk
    // Local horizontal edges are looped first
    for (int j=0; j<mi.MPIlocalVertexSize[1]; j++){
        for (int i=0; i<mi.MPIlocalCellSize[0]; i++){

//            computeEdgeFlux(mi, i, j, edgeflux, 0, mlu, gCellOut, gCellIn, 
//                            gwe, gpe, edgeEndsVertex, edgeEndsIndice,
//                            gauss_p, extractHoriEdgeInfo);

        }
    }

    // Loop vertical edges second
    for (int j=0; j<mi.MPIlocalCellSize[1]; j++){
        for (int i=0; i<mi.MPIlocalVertexSize[0]; i++) {

//            computeEdgeFlux(mi, i, j, edgeflux, mi.MPIlocalHoriEdgeSize, 
//                            mlu, gCellOut, gCellIn, 
//                            gwe, gpe, edgeEndsVertex, edgeEndsIndice,
//                            gauss_p, extractHoriEdgeInfo);

        }
    }

    return edgeflux;
}

double cellFlux(const MeshInfo& mi,
                const indice& lCell, 
                double * edgeFluxAdv,
                double * edgeFluxDif){

    // Return summation of flux on the edges of a given cell

    double flux = 0.0;

    // Four edge index and corresponding local index
    indice left   = lCell;
    indice right  = {lCell[0]+1, lCell[1]};
    indice bottom = lCell;
    indice top    = {lCell[0], lCell[1]+1};

    int left_flat  = FlatIndic(mi.MPIlocalVertexSize[0],left);
    int right_flat = FlatIndic(mi.MPIlocalVertexSize[0],left);
    int bottom_flat= FlatIndic(mi.MPIlocalCellSize[0], bottom) + mi.MPIlocalHoriEdgeSize;
    int top_flat   = FlatIndic(mi.MPIlocalCellSize[0], top) + mi.MPIlocalHoriEdgeSize;       

    flux = edgeFluxAdv[left_flat] - edgeFluxAdv[right_flat] + edgeFluxAdv[bottom_flat] - edgeFluxAdv[top_flat]; 

    return flux;
}
