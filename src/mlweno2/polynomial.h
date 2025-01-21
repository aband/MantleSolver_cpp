#ifndef POLYNOMIAL_H_
#define POLYNOMIAL_H_

#include "util.h"
#include "stencil.h"
#include <map>

// A 2D polynomial class
class polynomial {

    public:
        /*!
         * Highest polynomial order will be degree - 1
         */
        polynomial(const int& degreex, 
                   const int& degreey);

        ~polynomial();

        int setCoef (double * setcoef);

        /*!
         * Return a pointer pointing to a 
         * copy of polynomial coefficient.
         */
        double* getCoefPtr() {return coef;};

        /*!
         * Print out coefficients.
         */
        int printCoef() const;

        /*!
         * Reset polynomial degrees, which should only be used in testing
         */
        int resetDegree(const int& degreex,
                        const int& degreey);

        /*! 
         * Evaluation of point value with Horner's method
         */
        double eval(const double& x, const double& y) const;
        double operator() (const double& x, const double& y) const {return eval(x,y);}

         /*! 
         * Evaluation of point value with Horner's method
         */
        int evalDer(const int& derx,    const int& dery,
                    const int& degreex, const int& degreey,
                    const double& x,    const double& y);

    private:

        int degree[2] = {-1,-1};

        // Store coefficient in a 1D array
        double * coef = nullptr;
};

#endif
