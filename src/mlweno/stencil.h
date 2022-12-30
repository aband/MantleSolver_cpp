#ifndef STENCIL_H_
#define STENCIL_H_

/*
 *Define multi-level stencils for weno reconstruction
 *Designed for cartesian logically rectangular mesh.
 */

#include "util.h"

namespace MLWENO {

    enum stencilType {cellCentered, edgeCentered, vertexCentered};

    template <class T>
    class stencil{
        public:
            stencil(size_t I):                     I_(I), stencilSize_(I) {};
            stencil(size_t I, size_t J):           I_(I), J_(J), stencilSize_(I*J) {};
            stencil(size_t I, size_t J, size_t K): I_(I), J_(J), K_(K), stencilSize_(I*J*K) {};

            void CreateStencil() {stencil_.resize(stencilSize_);}

            T &operator()(size_t i){
                return stencil_[i];
            }

            T &operator()(size_t i, size_t j){
                return stencil_[i+j*I_];
            }

            T &operator()(size_t i, size_t j, size_t k){
                return stencil_[i+j*I_+k*I_*J_];
            }

            const vector<T> GetStencil() {return stencil_;}

            ~stencil();

        private:
            size_t I_;
            size_t J_;
            size_t K_; 
            size_t stencilSize_;
            vector<T> stencil_;
    };

// End of using name space MLWENO
} 

#endif
