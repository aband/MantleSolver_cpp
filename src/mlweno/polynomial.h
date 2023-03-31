#ifndef POLYNOMIAL_H_
#define POLYNOMIAL_H_

#include "util.h"
#include "stencil.h"

namespace tensorProductPoly{
    // Tensor product 2D polynomial for stencil basis
    class basePolynomial {

        public:
            basePolynomial() {};
            basePolynomial(const int maxDegree[2]);
            basePolynomial(const int maxDegree[2], double* coef);
            ~basePolynomial();
 
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
            stencilPolynomial(const indice& start, const vertex& center);
            stencilPolynomial(const indice& start, const vertex& center, 
                              const vector<indice>& targetCell);
            ~stencilPolynomial() {stencilPolyn_.clearPtr();
                                  delete collapsePolyn_;};

            void SetStencilPolynomials(const MeshInfo& mi,
                                       const stencil <indice>& stencilIndice);

            void SetUpScale(const MeshInfo& mi, const vector<indice>& targetCell) 
                           {SetTargetCell_(targetCell);ComputeCellBasedScale_(mi);};
            void SetUpScale(const MeshInfo& mi, const stencil <indice>& stencilIndice) {ComputeStencilBasedScale_(mi, stencilIndice);};

            const double GetScale() const {return scale_;}; 

            double eval(const double x, const double y) const;
            double eval(const vertex& P) const {return eval(P[0],P[1]);};
            double operator() (const double x, const double y) const {return eval(x,y);};
            double operator() (const vertex& P) const {return eval(P);};

            // Evaluation of point value for a given polynomial.
            // Separately without using collapsed polynomial.
            double eval(double x, double y, int poly) const;
            double eval(vertex P, int poly) const {return eval(P[0],P[1],poly);};
            double operator() (double x, double y, int poly) const {return eval(x,y,poly);};
            double operator() (vertex P, int poly) const {return eval(P,poly);};

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

            void SetTargetCell_(const vector<indice>& targetCell);
            vector<indice> targetCell_;
            double scale_ = -1.0;

            double coef_ = 1.0;

            // Calculate scale based on the information of a single cell
            void ComputeCellBasedScale_(const MeshInfo& mi);
            // Calculate scale based on the value of a given stencil
            void ComputeStencilBasedScale_(const MeshInfo& mi, const stencil <indice>& stencilIndice);

            stencil <basePolynomial*> stencilPolyn_;

            void SetCollapsePolyn_(const MeshInfo& mi, const stencil <indice>& stencilIndice);
            basePolynomial* collapsePolyn_ = nullptr;

            // Create polynomial smoothness indicator
            void CreateXi_();
            vector<double> Xi_;

            int maxR_;
            int minR_;

            void EvalSmoothIndic_(const MeshInfo& mi, const stencil <indice>& stencilIndice);
            double smoothnessIndic_ = -1;

    };

// End of using name space MLWENO
}

#endif
