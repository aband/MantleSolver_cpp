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

    if (OutBndry(mi, gCellIn)){
        return gCellOut;
    } else if (OutBndry(mi, gCellOut)){
        return gCellIn;
    } else {
        // Both gCellIn and gCellOut are inside boundary
        return gCellIn;
    }

}

inline double edgeFlux(const MeshInfo& mi,
                       const MLWENO::MLWENOUse& mluIn,
                       const MLWENO::MLWENOUse& mluOut,
                       const indice& globalCellIn,
                       const indice& globalCellOut,
                       const std::string& locationIn,
                       const std::string& locationOut,
                       const std::vector<vertex>& edge,
                       const std::vector<vertex>& gauss_p,
                       const std::vector<vertex>& velocity,
                       const std::valarray<double>& gwe,
                       const vetor<double>& param,
                       fluxFunc       fluxfuncAdv,
                       fluxFuncBndry  fluxfuncbndryAdv,
                       fluxFunc       fluxfuncDif,
                       fluxFuncBndry  fluxfuncbndryDif){
    // Compute integral of flux on edge
    // mluIn used for value where unit normal points to
    // mluOut used for value where unit normal tails at

    // Calculate unit normal vector and length of corresponding edge
    double len = length(edge);
    vertex unitNormal = UnitNormal(edge, len); 

    double flux = 0.0;
    vector<double> fluxpoint;
    vector<double> direction;

    for (const auto& it : velocity){
        direction.push_back(it[0]*unitNormal[0] + 
                            it[1]*unitNormal[1]); 
    }

    // Distinguish different boundary conditions
    if(OutBndryCell(mi, globalCellIn)){
        fluxpoint = fluxfuncbndryAdv(mi, mlu, edge, unitNormal, len, globalCellOut,
                                     locationOut, param, AssignBoundary(globalCellOut)) +
                    fluxfuncbndryDif(mi, mlu, edge, unitNormal, len, globalCellOut,
                                     locationOut, param, AssignBoundary(globalCellOut)) ;

    }else if (OutBndryCell(mi, globalCellIn)){
        fluxpoint = fluxfuncbndryAdv(mi, mlu, edge, unitNormal, len, globalCellIn,
                                     locationIn, param, AssignBoundary(globalCellIn)) + 
                    fluxfuncbndryDif(mi, mlu, edge, unitNormal, len, globalCellIn,
                                     locationIn, param, AssignBoundary(globalCellIn)) ;

    } else {
        // Interior edge
        fluxpoint = fluxfuncAdv(mi, mluIn, mluOut, edge, unitNormal, 
                                len, globalCellIn, globalCellOut,
                                locationIn, locationOut, gauss_p, gpe, param) + 
                    fluxfuncDif(mi, mluIn, mluOut, edge, unitNormal, 
                                len, globalCellIn, globalCellOut,
                                locationIn, locationOut, gauss_p, gpe, param) ;
    }

    // Integrate with gauess quadratrue scheme
    for (int g=0; g<gwe.size(); g++){
        flux += gwe[g] * fluxpoint.at(g) * direction.at(g); 
    }

    return flux;
}

double * edgeFluxAll(const MeshInfo* mi,
                     Vec * sol_darcy, Vec * g_darcy,
                     Vec * sol_stokes, Vec * g_stokes,
                     const int* refmap_darcy,
                     const int* refmap_stokes,
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

    vector<vertex> velocity_darcy;
    vector<vertex> velocity_stokes;
    vector<vertex> velocity_effect;

    // Interpolation of velocity
    velocity_darcy = ExtractVelocity(sol_darcy, g_darcy, refmap); 

    // Transform scaled variable to unscaled variable

    // Loop through the entire local mesh chunk
    // Local horizontal edges are looped first
    for (int j=0; j<mi.MPIlocalVertexSize[1]; j++){
        for (int i=0; i<mi.MPIlocalCellSize[0]; i++){

            // Flaten to integer
            indice local {i,j};
            int flatlocal = FlatIndic(mi.MPIlocalCellSize[0], local);

            // Identify location with global indices 
            indice gCellOut = local + mi.MPIlocalCellStart;
            indice gCellIn  = gCellOut - {0,1};

            // Extract edge vertex from mesh
            // Edge vertex from left to right
            start = local + ghostShift;
            end   = start + {1,0};

            start = mi.lmesh[FlatIndic(mi.MPIlocalVertexFull[0], start)]; 
            end   = mi.lmesh[FlatIndic(mi.MPIlocalVertexFull[0], end)];

            vector<vertex> edge {start, end};

            double len = length(edge);
            vertex unitNormal = UnitNormal(edge, len);

            // Get gauss quadrature points
            for (int g=0; g<gpe.size(); g++){gauss_p.at(g) = GaussMapPointsEdge({gpe[g]}, edge)};

            // Pick cell inside computational domain
            indice gcell_inside = PickCellInside(mi, gCellIn, gCellOut);

            // Compute velocity from computating results of stokes and darcy problems
            velocity_darcy = ExtractVelocity(sol_darcy, g_darcy, refmap_darcy, 
            mi, gauss_p, gcell_inside, hdiv, mybasis);

            velocity_stokes = ExtractVelocity(sol_stokes, g_stokes, refmap_stokes, 
            mi, gauss_p, gcell_inside, br, mybasis);

            // Compute effective velocity 

            edgeflux[flatlocal] = edgeFlux(mi, mluIn, mluOut, gCellIn, gCellOut,
            location(gCellIn), location(gCellOut), edge, gauss_p, velocity, gwe, {0.0},fluxfunc, fluxfuncbndry);

        }
    }

    // Loop vertical edges second
    for (int j=0; j<mi.MPIlocalCellSize[1]; j++){
        for (int i=0; i<mi.MPIlocalVertexSize[0]; i++) {
        
            indice local {i,j};
            int flatlocal = FlatIndic(mi.MPIlocalVertexSize[0], local);

            indice gCellOut = local + mi.MPIlocalCellStart;
            indice gCellIn  = gCellOut - {1,0};

            // Extract edge vertex from mesh
            // Edge vertex from top to bottom
            end = local + ghostShift; // bottom
            start = end + {0,1}; // top

            vector<vertex> edge {start, end};
            double len = length(edge);
            vertex unitNormal = UnitNormal(edge, len);

            // Get gauss quadrature points
            for (int g=0; g<gpe.size(); g++){gauss_p.at(g) = GaussMapPointsEdge({gpe[g]}, edge);}

            // Pick cell inside computational domain
            indice gcell_inside = PickCellInside(mi, gCellIn, gCellOut);

            // Compute velocity from computating results of stokes and darcy problems 
            velocity_darcy = ExtractVelocity(sol_darcy, g_darcy, refmap_darcy, 
            mi, gauss_p, gcell_inside, hdiv, mybasis);

            velocity_stokes = ExtractVelocity(sol_stokes, g_stokes, refmap_stokes, 
            mi, gauss_p, gcell_inside, br, mybasis);

            edgeflux[mi.MPIlocalHoriEdgeSize + flatlocal] = edgeFlux(mi, mluIn, mluOut, gCellIn, gCellOut,
            location(gCellIn), location(gCellOut), edge, gauss_p, velocity, gwe, {0.0},fluxfunc, fluxfuncbndry);

        }
    }

    return edgeflux;
}

double cellFlux(const indice& lCell, double * edgeFlux){

    // Return summation of flux on the edges of a given cell

    double flux = 0.0;

    // Four edge index and corresponding local index
    indice left   = lCell;
    indice right  = lCell + {1,0};
    indice bottom = lCell;
    indice top    = lCell + {0,1};

    int left_flat  = FlatIndic(mi.MPIlocalVertexSize[0],left);
    int right_flat = FlatIndic(mi.MPIlocalVertexSize[0],left);
    int bottom_flat= FlatIndic(mi.MPIlocalCellSize[0], bottom) + mi.localHoriEdgeSize;
    int top_flat   = FlatIndic(mi.MPIlocalCellSize[0], top) + mi.localHoriEdgeSize;       

    flux = edgeFlux[left_flat] + edgeFlux[right_flat] + edgeFlux[bottom_flat] + edgeFlux[top_flat]; 

    return flux;
}
