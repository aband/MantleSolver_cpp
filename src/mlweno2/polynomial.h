#ifndef POLYNOMIAL_H_
#define POLYNOMIAL_H_

#include "util.h"
#include "tensor.h"
#include <map>

// Some auxilliary funcitons defined for using polynomial class
int computeDerivative(const int& der,  const int& degree, 
                      const double& x, const double& scale,
                      double * coef, double * work);

// A 2D polynomial class
class polynomial {

    public:
        /*!
         * Highest polynomial order will be degree - 1
         */
        polynomial() {};

        polynomial(const int& degreex, 
                   const int& degreey);

        ~polynomial();

        int setCoef (double * setcoef);

        int setCoef (int index, double val);

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
        double eval(const vertex& p) const {return eval(p[0],p[1]);};
        double operator() (const vertex& p) const {return eval(p);}

         /*! 
         * Evaluation of point value with Horner's method
         */
        int evalDer(const int& derx,    const int& dery,
                    const double& x,    const double& y,
                    const double& scale, Tensor<double>& tensor);

    private:

        int degree[2] = {-1,-1};

        // Store coefficient in a 1D array
        double * coef = nullptr;
};

// Special numerical integral function
double polyNumIntegralFace(const vector<vertex>& corners,
                           const double& h,
                           const vertex& center,
                           polynomial& mypoly);

#endif
