#include "func.h"

/**
 * Change functions for transport part here.
 * Define transport functions and derivatives
 * 2D Burger's equation
 *
 */
double funcX(vertex x, double u, double t){
    return 0.5*u*u;
}

double dfuncX(vertex x, double u, double t){
    return u;
}

double funcY(vertex x, double u, double t){
    return 0.5*u*u;
}

double dfuncY(vertex x, double u, double t){
    return u;
}

/**
 * Define lax-friedrichs flux for transport problem.
 */

namespace LaxFriedrichs {

    double flux(const double& uIn, const double& uOut, const vertex& unitNormal, const vertex& point){
        double work = 0.0;

        work = (funcX(point, uIn, 0) + funcX(point, uOut, 0))*unitNormal[0] + 
               (funcY(point, uIn, 0) + funcY(point, uOut, 0))*unitNormal[1]; 

        /**
         * Using local lax friedrichs stabilization without passing global factor.
         */
        double alphaLF = max(fabs(dfuncX(point, uIn, 0)*unitNormal[0] + dfuncY(point, uIn, 0)*unitNormal[1]),
                             fabs(dfuncX(point, uOut,0)*unitNormal[0] + dfuncY(point, uOut,0)*unitNormal[1]));

        work = 0.5 * (work - alphaLF*(uOut - uIn));

        return work;
    }

    double flux(const double& uIn, const double& uOut, const vertex& unitNormal, const vertex& point, const double& alphaLF){
        double work = 0.0;
        
        work = (funcX(point, uIn, 0) + funcX(point, uOut, 0))*unitNormal[0] + 
               (funcY(point, uIn, 0) + funcY(point, uOut, 0))*unitNormal[1]; 

        /** 
         * Passing global lax friedrichs stabilization factor into this function.
         */

        work = 0.5 * (work - alphaLF*(uOut - uIn));

        return work; 
    }

    unordered_map<int,double> dflux(const double& uIn, const double& uOut, const vertex& unitNormal, 
                                    const vertex& mapped, const double& alphaLF, 
                                    const unordered_map<int, double>& duIn, 
                                    const unordered_map<int, double>& duOut){
        unordered_map<int, double> work;

        double in = 0.5*(dfuncX(mapped, uIn, 0)*unitNormal[0] + 
                         dfuncY(mapped, uIn, 0)*unitNormal[1] + alphaLF);

        double out = 0.5*(dfuncX(mapped, uOut, 0)*unitNormal[0] + 
                          dfuncY(mapped, uOut, 0)*unitNormal[1] - alphaLF);

        //! Loop through derivative of outside cell
        for (auto& duin: duIn){
            // No need to check if key exists for the fact that work is now completely empty
            work.insert(std::pair<int, double>(duin.first, duin.second*in));
        }

        //! Check if it is outside the boundary
        if (duOut.empty()==0) {
            //! Loop through derivative of inside cell
            for (auto& duout: duOut){
                if (work.count(duout.first) > 0){
                    // this key does exists
                    work[duout.first] += duout.second*out;
                } else {
                    // this key does not exists
                    work.insert(std::pair<int,double>(duout.first, duout.second*out));
                }
            }
        }

        return work;
    }

    double dflux(double uIn, double uOut, vertex unitNormal, vertex point){
        double work = 0.0;

        return work;
    }

}
