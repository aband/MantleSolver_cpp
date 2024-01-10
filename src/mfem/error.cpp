#include "error.h"

std::vector<double> GetFullSol(Vec * u, const bndryVal& bndryval, int M, int N){

    Vec sol = *u;

    double *arrayu;

    VecGetArray(sol, &arrayu);

    std::vector<double> work;

    work.resize(M*N);

    int count = 0;
    // Combine computed solution and restricted boundary values
    for (int j=0; j<N; j++){
        for (int i=0; i<M; i++){

            auto itFind = bndryval.find(j*M+i);
            if (itFind == bndryval.end()){
                work[j*M+i] = arrayu[count]; 
                count ++;
            } else {
                work[j*M+i] = itFind.second.DirichletVal;
            }
        }
    }

    VecRestoreArray(sol, &arrayu);

    return work;
}

std::array<double, 8> ExtractWeights(const std::vector<double>& fullsol, 
                                     const Hdivmixed& hdiv_,  
                                     const indice& globalElem,
                                     const MeshInfo& mi){

    std::array<double, 8> work;

    int flatGlobal = FlatIndic(mi, globalElem);

    std::array<int, 8> tmp = hdiv_.LocalToGlobal(mi, globalElem);

    for (int g=0; g<8; g++){
        work[g] = fullsol.at(tmp[g]);
    }

    return work;
}

double L2ErrorElemInterior(const std::array<double,8>& weight, 
                           const indice& globalElemIndic,
                           std::array<double,3> (*func)(const vertex& point),
                           const valarray<double>& gwf,
                           const vector<vertex>& gpf,
                           basis& basis_,
                           Hdivmixed& hdiv_){

    double elemError = 0.0;
    // Calculate the L2 Error on the given interior element 

    for (int g=0; g<gwf.size(); g++){
        //Loop through gauess quadrature points
        vertex mapped = GaussMapPointsFace(gpf[g], basis_.corners());
        double jac = abs(GaussJacobian(gpf[g], basis_.corners()));
        double gw = gwf[g];
        // Evaluate all eight basis functions for the target element
        std::array<vertex, 8> hdivwork = hdiv_.ComputeHdivmixed(basis_, mapped);
        // Combine these values with weights (calculated solution)
        valarray<double> approxVal = {0.0,0.0};
        for (int i=0; i<8; i++){
            approxVal += weight[i]*hdivwork[i]; 
        }
        // Get exact values
        std::array<double,3> trueSol = func(mapped);

        valarray<double> diff {0.0,0.0};

        diff[0] = approxVal[0] - trueSol[0];
        diff[1] = approxVal[0] - trueSol[0];

        elemError += gw*jac*(diff[0]*diff[0] + diff[1]*diff[1]);
    }

    return pow(elemError,0.5);
}
