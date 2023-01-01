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
            stencilPolynomial(const indice& start, const vertex& center, 
                              const vector<indice>& targetCell, const MeshInfo& mi);
            ~stencilPolynomial() {};

            void SetStencilPolynomials(const MeshInfo& mi,
                                       const vector<stencil <indice>>& stencilIndice);

        private:

            vertex center_;
            indice start_;
            vector<indice> targetCell_;
            double scale_ = -1.0;

            vector<stencil <basisPolynomial*>> stencilPolyn_;
    };

// End of using name space MLWENO
}

#endif
