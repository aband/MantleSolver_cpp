#ifndef POLYNOMIAL_H_
#define POLYNOMIAL_H_

#include "util.h"
#include "tensor.h"
#include <map>

// Some auxilliary funcitons defined for using polynomial class
int computeDerivative(const int& der,  const int& degree, 
                      const double& x, const double& scale,
                      double * coef, double * work);

class polynomial {

    public:
        /*!
         * Highest polynomial order will be degree - 1
         */
        polynomial() {};

        polynomial(const int& degreex, 
                   const int& degreey);

        ~polynomial() {};

        int setCoef (double * setcoef);

        int setCoef (int index, double val);

        double getCoef (int index) const {return coef[index];};

        /*!
         * Print out coefficients.
         */
        int printCoef() const;

        /*!
         * Reset polynomial degrees, which should only be used in testing
         * Be careful with it, number of coefficient may change due to change of 
         * degree of x and y.
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
                    const double& scale, Tensor<double>& tensor) const;

        int evalDer(const int& alpha, const double& x, const double& y,
                    const double& scale, Tensor<double>& der) const;

    private:

        vector<int> degree {-1,-1};

        // Store coefficient in a 1D array
        vector<double> coef;
};

// Special numerical integral function
double polyNumIntegralFace(const vector<vertex>& corners,
                           const double& h,
                           const vertex& center,
                           polynomial& mypoly);



#endif
