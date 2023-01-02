#ifndef RECONSTRUCTION_H_
#define RECONSTRUCTION_H_

/*
 *A ML weno reconstruction contains several information:
 *1. Center point
 *2. Stencils
 *3. max polynomial order
 *
 *What happens in the class of Reconstruction:
 *1. Compute smoothenss indicator
 *2. Compute non linear weights
 *3. Compute reconstruction values
 *4. Compute derivatives of reconstruction values
 */

#include <vector>
#include <pair>
#include <set>
#include <unordered_set>
#include <array>
#include <valarray>
#include <algorithm>
#include <numeric>
#include <memory>

#include "lapacke.h"
#include "integral.h"
#include <assert.h>


namespace MLWENO {

    class Reconstruction {

        public:
            Reconstruction(){};
            ~Reconstruction(){};

            void AddStencil(int[2] stencilSize, int[2] shift) {stencilSize_.push_back(stencilSize);
                                                               shift_.push_hack(shift);};

            void UpdateStencilSizeMax(int[2] max) {stencilSizeMax_[0] = max[0]; 
                                                   stencilSizeMax_[1] = max[1];};

        private:
            vector<int[2]> stencilSize_;
            int stencilSizeMax_[2] {-1,-1};
            vector<int[2]> shift_;

            vector<stencil <indice>> stenilIndice_;
            vector<stencilPolynomial> stncilPolyn_;
    }

}
#endif
