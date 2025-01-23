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
        int preparesigma(const vector<double>& area,
                         const vector<vector<vertex>>& cornerSet, 
                         const vertex& center, const double& scale);

    private:

        vector<int> size {-1,-1};

        Tensor<polynomial> tensorpoly;

        Tensor<double> sigma;

        // Calculate cell wise sigma, being called by function preparesigma
        double cellsigma(const double& area, const vector<vertex>& corners, 
                         const vertex& center, const double& scale, const int& index);

};

#endif
