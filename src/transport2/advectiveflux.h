#ifndef ADVECTIVEFLUX_H_
#define ADVECTIVEFLUX_H_

#include "mluse.h"
#include "trans_param.h"

// Numerical flux
inline double flux(const double& uin,  const double& uout, 
                   const double& fin, const double& fout,
                   const double& LF){

    return 0.5*(fout + fin - LF*(uout - uin));
}

double fluxintegral(const vertex& unitnormal,
                    const double& len,
                    const vector<double>& uIn,
                    const vector<double>& uOut,
                    const vector<vertex>& param);

double fluxintegralbndry();

#endif
