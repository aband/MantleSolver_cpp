#include "transport.h"
#include "transfunc.h" 

double LFflux(double fneg, double fpos, double uneg, double upos, double alpha){
    return 0.5*(fneg + fpos - alpha*(upos-uneg)); 
}

// advective flux computed at the interior edges
double TransportVariable::advflux_edge(const MeshInfo& mi, 
                                       const vector<vertex>& vel, 
                                       const vector<double>& uneg,
                                       const vector<double>& upos,
                                       const vector<double>& fneg,
                                       const vector<double>& fpos,
                                       const vertexSet& edge,
                                       bool localLF, double gLF){

    double work = 0.0;

    // Compute flux integration along the edge
    const valarray<double>& gwe = GaussWeightsEdge;
	 const valarray<double>& gpe = GaussPointsEdge;

    double len = length(edge);
	 vertex unitNormal = UnitNormal(edge, len);

    for (int g=0; g<gpe.size(); g++){
        vertex mapped = GaussMapPointsEdge({gpe[g]}, edge);

        if(localLF){
            gLF = abs(vel.at(g)[0]*unitNormal[0] + vel.at(g)[1]*unitNormal[1]);
            gLF = find_max(abs(dfdu(uneg.at(g))),abs(dfdu(upos.at(g))))*gLF;
        }

        work += gwe[g] * LFflux(fpos.at(g), fneg.at(g), upos.at(g), uneg.at(g), gLF) * len/2.0;
    }

    return work;
}

// advective flux computed at the boundary edges
double TransportVariable::advflux_edge(const MeshInfo& mi, 
                                       const vector<vertex>& vel, 
                                       const vector<double>& u,
                                       const vector<double>& f,
                                       const vertexSet& edge,
                                       bool localLF, double gLF){

    double work = 0.0;

    // Compute flux integration along the edge
    const valarray<double>& gwe = GaussWeightsEdge;
	 const valarray<double>& gpe = GaussPointsEdge;

    double len = length(edge);
	 vertex unitNormal = UnitNormal(edge, len);

    for (int g=0; g<gpe.size(); g++){
        vertex mapped = GaussMapPointsEdge({gpe[g]}, edge);

        if(localLF){
            gLF = abs(vel.at(g)[0]*unitNormal[0] + vel.at(g)[1]*unitNormal[1]);
            gLF = find_max(abs(dfdu(u.at(g))),abs(dfdu(u.at(g))))*gLF;
        }

        work += gwe[g] * LFflux(f.at(g), f.at(g), u.at(g), u.at(g), gLF) * len/2.0;
    }

    return work;
}
