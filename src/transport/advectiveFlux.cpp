#include "advectiveFlux.h"

inline double LFFlux(const vector<double>& u, 
                     const vector<double>& fu,
                     const double& LF){
    // fu[0] corresponding to fuIn , fu[1] corresponding to fuOut
    // u[0] corresponding to uIn , u[1] corresponding to uOut

    return 0.5*(fu[0]+fu[1] - LF*(u[1] - u[0]));
}

inline double LFFlux(const double& uin,  const double& uout, 
                     const double& fuin, const double& fuout, 
                     const double& LF){

    return 0.5*(fuout + fuin - LF*(uout - uin));
}

// Linear transport 
double advFlux(const valarray<double>& gwe,
               const vector<vertex>& vel, 
               const vector<double>& uIn, 
               const vector<double>& uOut,
               const vertex& unitnormal,
               const double& len){

    // LF flux

    double work = 0.0;

    for (int g=0; g<gwe.size(); g++){
        double f = vel.at(g)[0] * unitnormal[0] + vel.at(g)[1]*unitnormal[1];

        work += gwe[g] * len/2.0 * LFFlux(uIn.at(g), uOut.at(g), 
                                          uIn.at(g)*f, uOut.at(g)*f, 
                                          f);
    }

    return work;
}

double advFluxBndry(const MeshInfo& mi,
                    const valarray<double>& gwe,
                    const vector<vertex>& vel,
                    const vector<double>& u,
                    const vertex& unitnormal,
                    const double& len,
                    const indice& gCell,
                    const int& edgeflag,
                    const std::string& field){

    double work = 0.0;

    bndryTypeTrans bt = AssignBoundary(mi,gCell,edgeflag,field);

    vector<double> bu;
    bu.resize(u.size());

    switch (bt){
  
        case flux :
            work = bndryFluxAdv(gCell, edgeflag);
        break;

        case freeFlow:
            work = advFlux(gwe, vel, u, u, unitnormal, len); 
        break;

        case noFlow:
            work = 0.0;
        break;

        case dirichletTrans:
            if (field == "HD"){
                for (auto& it: bu){
                    it = 0.308114;
                }
                work = advFlux(gwe, vel, bu, bu, unitnormal, len);
            }else {
                for (auto& it: bu){
                     it = 0.1;
                 }
                work = advFlux(gwe, vel, bu, bu, unitnormal, len);
            }
        break;

        default:
            std::cout << "Boundary condition not defined properly ." << std::endl;
        break;
    }

//    cout << work << endl;;
    return work;
}

// =========== Implicit =================================
