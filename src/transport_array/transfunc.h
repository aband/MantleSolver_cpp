#ifndef TRANSFUNC_H_
#define TRANSFUNC_H_

#include "util.h"
#include "reconstruction.h"
#include "advectiveflux.h"
#include "diffusiveflux.h"

// This is the only costumized function in simulation
// Transport functions
double advfunc(const double& u, const vertex& vel, const vertex& unitnormal);

double dfdu(const double& u);

int dadvfunc(const derivative& du, const double& u, const vertex& vel, const vertex& unitnormal, derivative& work);

double inflow(const vertex& point, const vector<double>& param);

// Advective flux
int computeEdgeFlux(vector<double>& edgeflux, double t,
                    const MeshInfo& mi, double ** localvals,
						  const vector<reconstruction>& my_recon,
						  const vector<tensorstencilpoly>& sten_lg,
						  const vector<tensorstencilpoly>& sten_sm);

// Take into a velocity field that transport the phenomenon
int computeEdgeFlux(const vector<vertex>& velocityField,
                    vector<double>& edgeflux, double t,
                    const MeshInfo& mi, double ** localvals,
                    const vector<reconstruction>& my_recon,
                    const vector<tensorstencilpoly>& sten_lg,
                    const vector<tensorstencilpoly>& sten_sm);

#endif
