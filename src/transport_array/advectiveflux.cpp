#include "advectiveflux.h"

int advflux_all(const vector<reconstruction>& my_recon,
                const vector<tensorstencilpoly>& sten_lg,
                const vector<tensorstencilpoly>& sten_sm,
                double ** localvals,
					 vector<double>& advflux,
					 const MeshInfo& mi,
					 bool localLF,
					 double globalLF){

    int M = mi.MPIglobalCellSize[0];
    int N = mi.MPIglobalCellSize[1];

    // Clear flux vector
    advflux.clear();
    advflux.resize(M*(N+1) + N*(M+1));

    const valarray<double>& gwe = GaussWeightsEdge;
    const valarray<double>& gpe = GaussPointsEdge;

    vector<double> quadwts;
    vector<vertex> quadpts;

    double len = 0.0;
    vertex unitNormal;

    double work = 0.0;

    double LF = globalLF;

    // Interior horizonal edges 
    for (int j=0; j<N-1; j++){
    for (int i=0; i<M  ; i++){

        indice gcell {i, j};
        vertexSet corners = extractCorners(mi, gcell);

        int neg =     j*M + i;
        int pos = (j+1)*M + i;

        work = 0.0;

        // Get quad points
        vertexSet hori {corners.at(3), corners.at(2)};

        len = length(hori);
        unitNormal = UnitNormal(hori, len);

        for (int g=0; g<gpe.size(); g++){
            vertex mapped = GaussMapPointsEdge({gpe[g]}, hori);

            double uneg = my_recon.at(neg).eval(locvals, mapped, sten_lg, sten_sm);
            double upos = my_recon.at(pos).eval(locvals, mapped, sten_lg, sten_sm);

            vertex fneg = advfunc(uneg, vel.at(g), unitNormal);
            vertex fpos = advfunc(upos, vel.at(g), unitNormal);

            if (localLF) {
                LF = ads(vel.at(g)[0] * unitNormal[0] + vel.at(g)[1] * unitNormal[1]);
					 LF = find_max(abs(dadvfunc(uneg)), abs(dadvfunc(upos))) * LF;
            }

            work += gwe[g] * LFflux(fneg, fpos, uneg, upos, LF)* len/2.0; 
        }

        advflux.at(pos) = work;
    }}

    // On the boundary hori edge
    for (int i=0; i<M; i++){
        // Two edges
        vertexSet top {};
        indice top {};
        vertexSet corners = extractCorners(mi, top);

        indice top {};
        vertexSet corners = extractCorners(mi, bottom);

      
    }

    return 1;
}
