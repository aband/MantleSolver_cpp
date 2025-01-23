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
         * Obtain cell indice where counting starts.
         */
        int getStart(const indice& input){start = input;  return 1;};
        int getStart(const int& inputx,
                     const int& inputy){start[0] = inputx; start[1] = inputy; return 1;};

        /**!
         * Compute stencil polynomial coefficients
         */
        int setCoef(const MeshInfo& mi,
                    const Tensor<indice>& stencilindice,
                    const vertex& center, const double& scale);

    private:

        vector<int> size {-1,-1};

        Tensor<polynomial> tensorpoly;

        indice start {-1,-1};

};

#endif
