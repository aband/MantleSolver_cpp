#include "advectiveflux.h"

int testlink(){

    cout << "Can be linked" << endl;

    return 1;
}

static bool isoutflow(vertex vel, vertex normal){

    if(vel[0]*normal[0] + vel[1]*normal[1] > 0){
        return true;
    } else {
        return false;
    }
}

// Interior edge
double advflux_edge(const reconstruction& recon_neg,
					     const reconstruction& recon_pos,
                    const vector<tensorstencilpoly>& sten_lg,
                    const vector<tensorstencilpoly>& sten_sm,
                    double ** localvals,
						  const vertexSet& edge,
						  const vector<vertex>& vel,
						  bool localLF,
						  double gLF){

    double work = 0.0;

    // Compute flux integration along the edge
    const valarray<double>& gwe = GaussWeightsEdge;
	 const valarray<double>& gpe = GaussPointsEdge;

    double len = length(edge);
	 vertex unitNormal = UnitNormal(edge, len);

    for (int g=0; g<gpe.size(); g++){
        vertex mapped = GaussMapPointsEdge({gpe[g]}, edge);

        double uneg = recon_neg.eval(localvals, mapped, sten_lg, sten_sm); 
        double upos = recon_neg.eval(localvals, mapped, sten_lg, sten_sm); 

        double fneg = advfunc(uneg, vel.at(g), unitNormal);
        double fpos = advfunc(upos, vel.at(g), unitNormal);

        if (localLF) {
            gLF = abs(vel.at(g)[0]*unitNormal[0] + vel.at(g)[1]*unitNormal[1]);
				gLF = find_max(abs(dfdu(uneg)),abs(dfdu(upos)))*gLF;
        }
        work += gwe[g] * LFflux(fneg, fpos, uneg, upos, gLF)* len/2.0;
    }

    return work;
}

// Boundary flux integration, use one-sided reconstruction object
// Used in free outflow boundary condition
double advflux_edge(const reconstruction& recon,
                    const vector<tensorstencilpoly>& sten_lg,
                    const vector<tensorstencilpoly>& sten_sm,
                    double ** localvals,
						  const vertexSet& edge,
						  const vector<vertex>& vel,
						  bool localLF,
						  double gLF){

    double work = 0.0;

    // Compute flux integration along the edge
    const valarray<double>& gwe = GaussWeightsEdge;
	 const valarray<double>& gpe = GaussPointsEdge;

    double len = length(edge);
	 vertex unitNormal = UnitNormal(edge, len);

    for (int g=0; g<gpe.size(); g++){
        vertex mapped = GaussMapPointsEdge({gpe[g]}, edge);

        double u = recon.eval(localvals, mapped, sten_lg, sten_sm); 
        double f = advfunc(u, vel.at(g), unitNormal);
 
        if (localLF) {
            gLF = abs(vel.at(g)[0]*unitNormal[0] + vel.at(g)[1]*unitNormal[1]);
				gLF = find_max(abs(dfdu(u)),abs(dfdu(u)))*gLF;
        }
        work += gwe[g] * LFflux(f, f, u, u, gLF)* len/2.0;
 
    }

    return work;
}

double advflux_edge(double (*func)(const vertex& point,
								           const vector<double>& param),
					     const vector<double>& param,
					     const vertexSet& edge,
						  const vector<vertex>& vel,
						  bool localLF,
						  double gLF){

    double work = 0.0;

    // Compute flux integration along the edge
    const valarray<double>& gwe = GaussWeightsEdge;
	 const valarray<double>& gpe = GaussPointsEdge;

    double len = length(edge);
	 vertex unitNormal = UnitNormal(edge, len);

    for (int g=0; g<gpe.size(); g++){
        vertex mapped = GaussMapPointsEdge({gpe[g]}, edge);

        double u = func(mapped, param);
        double f = advfunc(u, vel.at(g), unitNormal);
 
        if (localLF) {
            gLF = abs(vel.at(g)[0]*unitNormal[0] + vel.at(g)[1]*unitNormal[1]);
				gLF = find_max(abs(dfdu(u)),abs(dfdu(u)))*gLF;
        }
        work += gwe[g] * LFflux(f, f, u, u, gLF)* len/2.0;
 
    }

    return work;
}
