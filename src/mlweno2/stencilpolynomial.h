#ifndef STENCILPOLYNOMIAL_H_
#define STENCILPOLYNOMIAL_H_

#include "polynomial.h"

class stencilpolynomial {

    public:
        stencilpolynomial() {};

        stencilpolynomial(const int& sizex,
                          const int& sizey);

        ~stencilpolynomial() {};

        /**!
         * Compute stencil polynomial coefficients
         * Takes in one vector of all corners.
         * No direct involvement of MeshInfo object
         */
        int setCoef(const vector<vector<vertex>>& cornerSet,
                    const vertex& center, const double& scale);

        /**!
         * Print stencil polynomial coefficient
         */
        int printCoef();

        /**!
         * Evaluate with given points.
         */
        double eval(const Tensor<double>& sol, const vertex& point, const vertex& center, const double& h); 

        /**!
         * Evaluate smoothness indicator tensor.
         */
        int sigma(const vector<vertex>& corners, const double& area,
                  const vertex& center, const double& h);

    private:

        vector<int> size {-1,-1};

        Tensor<polynomial> tensorpoly;

        Tensor<double> tensorsigma;

        // Calculate the integral in smoothness indicator
        double sigmaintegral(const vector<vertex>& corners, const double& area,
                             const vertex& center, const double& h,
                             const int& index1, const int& index2);

};

#endif
