#include "advectiveflux.h"

double fluxintegral(const vertex& unitnormal,
                    const double& len,
                    const vector<double>& uin,
                    const vector<double>& uout,
                    const vector<vertex>& param){

    double work = 0.0;

    const valarray<double>& gwe = GaussWeightsEdge;

    vector<double> fin;  fin.resize(gwe.size());
    vector<double> fout; fout.resize(gwe.size());
    vector<double> LF;   LF.resize(gwe.size());

    advfunc(uin, uout, fin, fout, param, LF, unitnormal);

    for (int g=0; g<gwe.size(); g++){
        work += gwe[g] * len/2.0 * flux(uin.at(g), uout.at(g), fin.at(g), fout.at(g), LF.at(g)) * 
                (unitnormal[0]*param.at(0)[0] + unitnormal[1]*param.at(0)[1]);     
    }

    return work;
}

// For testing purpose right now
// Simple Reimann shock and rarefaction problem where no flux penetrates through boundary
double fluxintegralbndry(){

    return 0;
}
