#include "edgeFlux.h"

double edgeFlux(const MeshInfo& mi,
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
    vector<double> direction;

    for (const auto& it : velocity){
        direction.push_back(it[0]*unitNormal[0] + 
                            it[1]*unitNormal[1]); 
    }

    // Distinguish different boundary conditions
    if(OutBndry(globalCellIn)){
        fluxpoint = fluxfuncbndry(mi, mlu, edge, unitNormal, len, globalCellIn,
                                  locationIn, param, AssignBoundary(globalCellIn));

    }else if (OutBndry(globalCellOut)){
        fluxpoint = fluxfuncbndry(mi, mlu, edge, unitNormal, len, globalCellOut,
                                  locationOut, param, AssignBoundary(globalCellOut));
    } else {
        // Interior edge
        fluxpoint = fluxfunc(mi, mluIn, mluOut, edge, unitNormal, 
                             len, globalCellIn, globalCellOut,
                             locationIn, locationOut, gauss_p, gpe, param);
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
                     basis& mybasis,
                     BRMixed& br,
                     Hdivmixed& hdiv,
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

            for (int g=0; g<gpe.size(); g++){gauss_p.at(g) = GaussMapPointsEdge({gpe[g]}, edge)};

            // Compute velocity of from computation results of stokes and darcy problems
            velocity_darcy = ExtractVelocity(sol_darcy, g_darcy);
            velocity_stokes = ExtractVelocity(sol_stokes, g_stokes);

            // Compute effective velocity 

            edgeflux[flatlocal] = edgeFlux(mi, mluIn, mluOut, gCellIn, gCellOut,
            location(gCellIn), location(gCellOut), edge, gauss_p, velocity, gwe, {0.0},gluxfunc, fluxfuncbndry);

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

            for (int g=0; g<gpe.size(); g++){gauss_p.at(g) = GaussMapPointsEdge({gpe[g]}, edge);}


        }
    }

    return edgeflux;
}
