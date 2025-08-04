#ifndef STENCILPOLYNOMIAL_H_
#define STENCILPOLYNOMIAL_H_

#include "polynomial.h"

class stencilpolynomial {

    public:
        stencilpolynomial() {};

        stencilpolynomial(const int& sizex,
                          const int& sizey);

        stencilpolynomial(const int& order);

        ~stencilpolynomial() {};

        vertex center;
        double h;

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
        double eval(const Tensor<double>& sol, const vertex& point, const vertex& center, const double& h) const; 

        double eval(const Tensor<double>& sol, const vertex& point) const 
        {return eval(sol, point, center, h);};

        /*!
         * Evaluate base polynomial
         */
        double eval(const vertex& point, const vector<int>& baseindex) const;

        /**!
         * Evaluate smoothness indicator tensor.
         */
        int sigma(const vector<vertex>& corners, const double& area,
                  const vertex& center, const double& h);

        int sigmacomplete(const vector<vertex>& corners,
                          const double& area,
                          const vertex& center,
                          const double& h);

        double sigma(const Tensor<double>& sol) const;

        /**!
         * Compute derivative of smoothness indicator at the same time
         * Used when doing implicit time stepping
         */
        int sigma(const Tensor<double>& sol, double& sig, vector<double>& der);

        /**!
         * Get derivative of smoothness indicator only
         */
        int dsigma(const Tensor<double>& sol, vector<double>& der);

        /**!
         * Function will be used for testing purpose
         */
        polynomial createCollapsePoly(const Tensor<double>& sol);

        double sigma(const polynomial& collapse, const vector<vertex>& corners,
                     const double& area, const vertex& center, const double& h);

        int printSigmaTensor();

    private:

        vector<int> size {-1,-1};

        int order = 0;

        Tensor<polynomial> tensorpoly;

        Tensor<double> tensorsigma;

        // Calculate the integral in smoothness indicator
        double sigmaintegral(const vector<vertex>& corners, const double& area,
                             const vertex& center, const double& h,
                             const int& index1, const int& index2);

};

#endif
