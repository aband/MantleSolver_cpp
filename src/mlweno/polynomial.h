#ifndef POLYNOMIAL_H_
#define POLYNOMIAL_H_

#include "util.h"
#include "stencil.h"

namespace MLWENO{
    // Tensor product 2D polynomial for stencil basis
    class basisPolynomial {

        public:
            basisPolynomial() {};
            basisPolynomial(const int maxDegree[2]);
            basisPolynomial(const int maxDegree[2], double* coef);
            ~basisPolynomial();
 
            void setMaxDegree(const int maxDegree[2]);
            void setCoef(double* coef);

            // Evaluation of point value for a given polynomial
            double eval(double x, double y) const;
            double eval(vertex P) const {return eval(P[0],P[1]);};
            double operator() (double x, double y) const {return eval(x,y);};
            double operator() (vertex P) const {return eval(P);};

            // Evaluation of derivative value for a given polynomial
            double evalDer();

        private:

            int maxDegree_[2] = {-1,-1}; 
            double* coef_ = nullptr;
    };

    // stencil polynomials used for weno reconstruction
    class stencilPolynomial {

        public:
            stencilPolynomial() {};
            stencilPolynomial(const MeshInfo& mi, const vertex& center);
            ~stencilPolynomial() {};

            void SetStencilPolynomials(const vector<stencil <indice>>& stencilIndice);

        private:

            vertex center_;

            int stencilSize_[2];
            double scale_;

            stencil <basisPolynomial *> stencilPolyn_();
    };

// End of using name space MLWENO
}

#endif
