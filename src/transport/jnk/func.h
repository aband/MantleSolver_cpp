#ifndef FUNC_H_
#define FUNC_H_

#include "reconstruction.h"

//! Two dimension functions
//! Functions for different dimensions are realized via function overloading.
double funcX(vertex x, double u, double t);

double dfuncX(vertex x, double u, double t);

double funcY(vertex x, double u, double t);

double dfuncY(vertex x, double u, double t);

namespace LaxFriedrichs {
    //! Local lax friedrichs scheme
    double flux(const double& uIn, const double& uOut, 
                const vertex& unitNormal, const vertex& point);
    //! Global lax friedrichs scheme
    double flux(const double& uIn, const double& uOut, 
                const vertex& unitNormal, const vertex& point, const double& alphaLF);

    //! derivative of global lax friedrichs scheme
    unordered_map<int, double> dflux(const double& uIn, const double& uOut, const vertex& unitNormal, 
                                     const vertex& mapped, const double& alphaLF, 
                                     const unordered_map<int, double>& duOut, const unordered_map<int, double>& duIn);

}

#endif
