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

            double* getCoef() const;

            // Return the size of the polynomial
            int getSize() const {return maxDegree_[0]*maxDegree_[1];};

            // Return the size of degree x as I
            int getI() const {return maxDegree_[0];};
            // Return the size of degree y as J
            int getJ() const {return maxDegree_[1];};

            // Evaluation of point value for a given polynomial
            double eval(double x, double y) const;
            double eval(vertex P) const {return eval(P[0],P[1]);};
            double operator() (double x, double y) const {return eval(x,y);};
            double operator() (vertex P) const {return eval(P);};

            // Evaluation of derivative value for a given polynomial
            double der(int derX, int derY, double x, double y) const;
            double der(int derX, int derY, vertex P) const {return der(derX, derY, P[0], P[1]);};

            // Print coefficients out
            void printCoef() const; 

            double getCoef(int i) const;

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
            ~stencilPolynomial() {stencilPolyn_.clearPtr();
                                  delete [] collapsePolyn_;};

            void SetStencilPolynomials(const MeshInfo& mi,
                                       const stencil <indice>& stencilIndice);

            void SetUpScale() {ComputeCellBasedScal_();};
            void SetUpScale(const stencil <indice>& stencilIndice) {ComputeStencilBasedScale_(stencilIndice);};

            const double GetScale() const {return scale_;}; 

            double eval(const double x, const double y) const;
            double eval(const vertex& P) const {return eval(P[0],P[1]);};
            double operator() (const double x, const double y) const {return eval(x,y);};
            double operator() (const vertex& P) const {return eval(P);};

            // Check basis polynomial coefficients
            void printCoef();
            void printCoef(int s);

            // Set collapse polynomial
            void SetCollapsePolyn(const MeshInfo& mi, const stencil <indice>& stencilIndice);

            double GetSmoothIndic(const MeshInfo& mi, const stencil <indice>& stencilIndice);

            const int GetOrderX() const { return stencilPolyn_.getI();};
            const int GetOrderY() const { return stencilPolyn_.getJ();};

        private:

            vertex center_ = {0.0,0.0};
            indice start_ = {-1,-1};
            vector<indice> targetCell_;
            double scale_ = -1.0;

            // Calculate scale based on the information of a single cell
            void ComputeCellBasedScale_();
            // Calculate scale based on the value of a given stencil
            void ComputeStencilBasedScale_();

            stencil <basisPolynomial*> stencilPolyn_;

            void SetCollapsePolyn_(const MeshInfo& mi, const stencil <indice>& stencilIndice);
            basisPolynomial* collapsePolyn_ = nullptr;

            // Create polynomial smoothness indicator
            void EvalSmoothIndic_(const MeshInfo& mi, const stencil <indice>& stencilIndice);
            double smoothnessIndic_ = -1;
    };

// End of using name space MLWENO
}

#endif
