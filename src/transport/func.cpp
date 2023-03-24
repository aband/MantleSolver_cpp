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

    double flux(double uIn, double uOut, vertex unitNormal, vertex point){
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

    double flux(double uIn, double uOut, vertex unitNormal, vertex point, double alphaLF){
        double work = 0.0;
        
        work = (funcX(point, uIn, 0) + funcX(point, uOut, 0))*unitNormal[0] + 
               (funcY(point, uIn, 0) + funcY(point, uOut, 0))*unitNormal[1]; 

        /** 
         * Passing global lax friedrichs stabilization factor into this function.
         */

        work = 0.5 * (work - alphaLF*(uOut - uIn));

        return work; 
    }

    double dflux(double uIn, double uOut, vertex unitNormal, vertex point, double alphaLF){
        double work = 0.0;

        return work;
    }

    double dflux(double uIn, double uOut, vertex unitNormal, vertex point){
        double work = 0.0;

        return work;
    }

}
