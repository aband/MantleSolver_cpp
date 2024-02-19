#include "error.h"

std::vector<double> GetFullSol(Vec * u, const bndryVal& bndryval, int dof){

    Vec sol = *u;

    double *arrayu;

    VecGetArray(sol, &arrayu);

    std::vector<double> work;

    work.resize(floor(dof));

    int count = 0;
    // Combine computed solution and restricted boundary values
    for (int j=0; j<dof; j++){
            auto itFind = bndryval.find(j);
            if (itFind == bndryval.end()){
                work[j] = arrayu[count]; 
                count ++;
            } else {
                work[j] = itFind->second.DirichletVal;
            }
    }

    VecRestoreArray(sol, &arrayu);

    return work;
}

std::array<double, 8> ExtractWeights(const std::vector<double>& fullsol, 
                                     std::array<int, 8> ltgMap){

    std::array<double, 8> work;

    for (int g=0; g<8; g++){
        work[g] = fullsol.at(ltgMap[g]);
    }

    return work;
}

std::array<double, 12> ExtractWeights(const std::vector<double>& fullsol, 
                                      std::array<int, 12> ltgMap){
    std::array<double, 12> work;

    for (int g=0; g<12; g++){
        work[g] = fullsol.at(ltgMap[g]);
    }

    return work;
}

double L2ErrorElem(const std::array<double,8>& weight, 
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

//        cout << "Evaluation on gauss points : " << approxVal[0] << " " << approxVal[1] << endl;

        // Get exact values
        std::array<double,3> trueSol = func(mapped);

        valarray<double> diff {0.0,0.0};

        //diff[0] = approxVal[0] - trueSol[0];
        //diff[1] = approxVal[1] - trueSol[1];

        diff[0] = abs(approxVal[0] - trueSol[0]);
        diff[1] = abs(approxVal[1] - trueSol[1]);

        elemError += gw*jac*(diff[0]*diff[0] + diff[1]*diff[1]);
    }

    return elemError;
}

double L2ErrorElem(const std::array<double, 12>& weight,
                   const indice& globalElemIndic,
                   std::array<double, 3>(*func)(const vertex& point),
                   const valarray<double>& gwf,
                   const vector<vertex>& gpf,
                   basis& basis_,
                   BRMixed& br_){

    double elemError = 0.0;

    for (int g=0; g<gwf.size(); g++){
        // Loop through gauess quadrature points
        vertex mapped = GaussMapPointsFace(gpf[g], basis_.corners());
        double jac = abs(GaussJacobian(gpf[g], basis_.corners()));
        double gw = gwf[g];
        // Evaluate all eight basis functions for the target element
        std::array<vertex, 12> brwork = br_.ComputeBRmixed(basis_, mapped);
        // Combine these values with weights (calculated solution)
        valarray<double> approxVal = {0.0,0.0};
        for (int i=0; i<12; i++){
            approxVal += weight[i]*brwork[i];
        }

        // Get exact values
        std::array<double, 3> trueSol = func(mapped);
        valarray<double> diff {0.0,0.0};

        diff[0] = abs(approxVal[0] - trueSol[0]);
        diff[1] = abs(approxVal[1] - trueSol[1]);

        elemError += gw*jac*(diff[0]*diff[0] + diff[1]*diff[1]);
    }

    return elemError;
}

// Compute L2 error for pressure
double L2ErrorElem(const double& approxP,
                   std::array<double,3> (*func)(const vertex& point),
                   const valarray<double>& gwf,
                   const vector<vertex>& gpf,
                   basis& basis_,
                   const double& area){

    double elemError = 0.0;

    double elemSum = 0.0;

    for (int g=0; g<gwf.size(); g++){
        vertex mapped = GaussMapPointsFace(gpf[g], basis_.corners());

        double jac = abs(GaussJacobian(gpf[g], basis_.corners()));
        double gw = gwf[g];

        std::array<double, 3> trueP = func(mapped);

        elemError += gw*jac*pow(trueP[2] - approxP,2);
    }

    for (int g=0; g<gwf.size(); g++){
        vertex mapped = GaussMapPointsFace(gpf[g], basis_.corners());

        double jac = abs(GaussJacobian(gpf[g], basis_.corners()));
        double gw = gwf[g];

        std::array<double, 3> trueP = func(mapped);

        elemSum += gw*jac*trueP[2];
    }

    return elemError;
}

double L2ErrorElem(const std::array<double,8>& weight, 
                   const indice& globalElemIndic,
                   vertex (*func)(const vertex& point, PhysProperty * pp),
                   PhysProperty * pp,
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

//        cout << "Evaluation on gauss points : " << approxVal[0] << " " << approxVal[1] << endl;

        // Get exact values
        vertex trueSol = func(mapped,pp);

        valarray<double> diff {0.0,0.0};

        //diff[0] = approxVal[0] - trueSol[0];
        //diff[1] = approxVal[1] - trueSol[1];

        diff[0] = abs(approxVal[0] - trueSol[0]);
        diff[1] = abs(approxVal[1] - trueSol[1]);

        elemError += gw*jac*(diff[0]*diff[0] + diff[1]*diff[1]);
    }

    return elemError;
}

double L2ErrorElem(const std::array<double, 12>& weight,
                   const indice& globalElemIndic,
                   vertex (*func)(const vertex& point, PhysProperty * pp),
                   PhysProperty * pp,
                   const valarray<double>& gwf,
                   const vector<vertex>& gpf,
                   basis& basis_,
                   BRMixed& br_){

    double elemError = 0.0;

    for (int g=0; g<gwf.size(); g++){
        // Loop through gauess quadrature points
        vertex mapped = GaussMapPointsFace(gpf[g], basis_.corners());
        double jac = abs(GaussJacobian(gpf[g], basis_.corners()));
        double gw = gwf[g];
        // Evaluate all eight basis functions for the target element
        std::array<vertex, 12> brwork = br_.ComputeBRmixed(basis_, mapped);
        // Combine these values with weights (calculated solution)
        valarray<double> approxVal = {0.0,0.0};
        for (int i=0; i<12; i++){
            approxVal += weight[i]*brwork[i];
        }

        // Get exact values
        vertex trueSol = func(mapped,pp);
        valarray<double> diff {0.0,0.0};

        diff[0] = abs(approxVal[0] - trueSol[0]);
        diff[1] = abs(approxVal[1] - trueSol[1]);

        elemError += gw*jac*(diff[0]*diff[0] + diff[1]*diff[1]);
    }

    return elemError;
}


