#include "diffusiveFlux.h"
#include "lagrange_tmp.h"

/**!
 * Assign values evaluating diffusion functions along the normal direction.
 */
inline void assignDiffVals(const MLWENO::MLWENOUse& mlu,
                           const MeshInfo&mi,
                           const indice& globalCellIn,
                           const indice& globalCellOut,
                           const int& locationIn,
                           const int& locationOut,
                           const vertex& unitNormal,
                           const vertex& mapped,
                           const double& dx,
                           vector<double>& diffVals){

    int halfnumPts = diffVals.size() / 2;

    for (int i=0; i<halfnumPts; i++){
        vertex point0 = mapped - (halfnumPts - 0.5)*dx*unitNormal;
        vertex point1 = mapped + 0.5*dx*unitNormal;

        diffVals[i] = diffFunc(mlu.Evaluate(point0,globalCellIn,mi,locationIn));
        diffVals[i+halfnumPts] = diffFunc(mlu.Evaluate(point1,globalCellOut,mi,locationOut));
    }

}

double getDifFlux(const MLWENO::MLWENOUse& mlu,
                  const MeshInfo& mi,
                  const std::array<vertex,2>& edge,
                  const vertex& unitNormal,
                  const double& len,
                  const indice& globalCellIn,
                  const indice& globalCellOut,
                  const int& locationIn,
                  const int& locationOut,
                  const valarray<double>& gwe,
                  const valarray<double>& gpe,
                  const double& scale){

    double work = 0.0;

    int degree = gwe.size() + 1; 
    const int numPts = std::ceil((degree+1)/2.0) * 2;

    // Get diameter
    const double hIn  = mi.cellArea.at(FlatIndic(mi,globalCellIn));
    const double hOut = mi.cellArea.at(FlatIndic(mi,globalCellOut));
    
    const double h = scale * ((hIn < hOut) ? hIn : hOut);

    const double dx = h /(double)(numPts - 1);

    LagrangeBasisDeriv lagDer(numPts - 1);

    const int halfPts = numPts/2;

    // Diffusion function are evaluated at the points in the order of 
    // from outside to inside.
    vector<double> diffVals(numPts, 0);

    vertexSet tmpEdge = {edge[0], edge[1]};

    for (int g=0; g<gwe.size(); g++){
        assignDiffVals(mlu, mi, globalCellIn, globalCellOut, locationIn, 
                       locationOut, unitNormal, GaussMapPointsEdge({gpe[g]}, tmpEdge), 
                       dx, diffVals);

        for (int i=0; i<numPts; i++){
            work -= lagDer.middle(numPts-1, i) / dx * gwe[g] * len/2.0 *
                    diffVals.at(i);
        }
    }

    return work;
}

double getDifFlux(const MLWENO::MLWENOUse& mlu,
                  const MeshInfo& mi,
                  const std::array<vertex,2>& edge,
                  const vertex& unitNormal,
                  const double& len,
                  const indice& globalCellIn,
                  const int& locationIn,
                  const valarray<double>& gwe,
                  const valarray<double>& gpe,
                  const double& scale,
                  bndryTypeTrans bt){

    double work = 0.0;

    int degree = gwe.size() + 1;
    const int numPts = std::ceil((degree+1)/2.0) * 2;

    // No outside cell for cells on the boundary
    const double h = mi.cellArea.at(FlatIndic(mi, globalCellIn));

    const double dx = h /(double)(numPts - 1);

    LagrangeBasisDeriv lagDer(numPts - 1);

    const int halfPts = numPts/2;

    // Diffusion function are evaluated at the points in the order of 
    // from outside to inside.
    vector<double> diffVals(numPts, 0);

    vertexSet tmpEdge = {edge[0], edge[1]};



    return work;
}

// ================= new standard compatible diffusive flux function
// Being regarded as a fluxfunc type of function
double difFlux(const valarray<double>& gwe,
               const vector<vertex>& vel,
               const vector<double>& uIn,
               const vector<double>& uOut,
               const vertex& unitnormal,
               const double& len){

}

