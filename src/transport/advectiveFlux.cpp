#include "advectiveFlux.h"

inline double LFFlux(const vector<double>& u, 
                     const vector<double>& fu,
                     const double& LF){
    // fu[0] corresponding to fuIn , fu[1] corresponding to fuOut
    // u[0] corresponding to uIn , u[1] corresponding to uOut

    return 0.5*(fu[0]+fu[1] - LF*(u[1] - u[0]));
}

inline double LFFlux(const double& uin, const double& uout, const double& fuin, const double& fuout, const double& LF){

    return 0.5*(fuout + fuin - LF*(uout - uin));
}

vector<double> advFlux(const MeshInfo& mi,
                       const MLWENO::MLWENOUse& mlu,
                       const vector<vertex>& edge,
                       const vertex& unitNormal,
                       const double& len,
                       const indice& gCellIn,
                       const indice& gCellOut,
                       const std::string& locIn,
                       const std::string& locOut,
                       const vector<double>& LFparam,
                       const vector<double>& direction,
                       const valarray<double>& gpe){

    assert(LFparam.size() == gpe.size());

    vector<double> work;
    work.resize(gpe.size()); 

    for (int g=0; g<gpe.size(); g++){

        vertex mapped = GaussMapPointsEdge({gpe[g]}, edge);
        double uIn  = mlu.Evaluate(mapped, gCellIn, mi, locIn);
        double uOut = mlu.Evaluate(mapped, gCellOut, mi, locOut);

        // Compute Lax_Friedrich flux
        // Here transport equation is hard coded within
        work.at(g) = LFFlux({uIn, uOut}, 
          {direction.at(g)*uIn, direction.at(g)*uOut}, 
          LFparam.at(g));

    }

    return work;
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

        work += gwe[g] * len/2.0 * LFFlux(uIn.at(g), uOut.at(g), uIn.at(g)*f, uOut.at(g)*f, 
                                          f);
    }

    return work;
}

vector<double> advFluxBndry(const MeshInfo& mi,
                            const MLWENO::MLWENOUse& mlu,
                            const vector<vertex>& edge,
                            const vertex& unitNormal,
                            const double& len,
                            const indice& gCell,
                            const std::string& loc,
                            const vector<double>& direction,
                            const vector<double>& LFparam,
                            const valarray<double>& gpe,
                            const bndryTypeTrans& bt, 
                            const int& flag,
                            const int& locedge){

    assert(LFparam.size() == gpe.size());

    vector<double> work;
    work.resize(gpe.size()); 

    vector<double> fu;
    vector<double> u;

    fu.resize(2);
    u.resize(2);
   
    for (int g=0; g<gpe.size(); g++){

        vertex mapped = GaussMapPointsEdge({gpe[g]}, edge);

        switch (bt){

            case flux:

                work.at(g) = bndryFluxAdv(gCell, locedge);

            break;

            case dirichletTrans:

                u.at(flag)   = mlu.Evaluate(mapped, gCell, mi, loc);
                u.at(1-flag) = bndryValAdv(gCell, locedge);

                fu.at(flag)   = direction.at(g)*u.at(flag);
                fu.at(1-flag) = direction.at(g)*u.at(1-flag);

                work.at(g) = LFFlux(u, fu, LFparam.at(g));

            break;

            case freeFlow:

                u.at(flag) = mlu.Evaluate(mapped, gCell, mi, loc);
                u.at(1-flag) = u.at(flag);

                fu.at(flag)   = direction.at(g)*u.at(flag);
                fu.at(1-flag) = direction.at(g)*u.at(1-flag);

                work.at(g) = LFFlux(u, fu, LFparam.at(g));

            break;

            default:

                std::cout << "Boundary condition not defined properly ." << std::endl;

            break;
        }

    }

    return work;
}

double advFluxBndry(const valarray<double>& gwe,
                    const vector<vertex>& vel,
                    const vector<double>& u,
                    const vertex& unitnormal,
                    const double& len,
                    const bndryTypeTrans& bt,
                    const int& edgetype,
                    const indice& gCell){

    double work = 0.0;

    switch (bt){
  
        case flux :
            work = bndryFluxAdv(gCell, edgetype);
        break;

        case freeFlow:
            work = advFlux(gwe, vel, vel, u, u, unitnormal, len); 
        break;

        default:
            std::cout << "Boundary condition not defined properly ." << std::endl;
        break;
    }

    return work;
}

// =========== Implicit =================================
