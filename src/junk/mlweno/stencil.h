#ifndef STENCIL_H_
#define STENCIL_H_

/*
 *Define multi-level stencils for weno reconstruction
 *Designed for cartesian logically rectangular mesh.
 */

#include "util.h"

enum stencilType {cellCentered, edgeCentered, vertexCentered};

template <class T>
class stencil{
    public:
        stencil() {};
        stencil(size_t I):                     I_(I), stencilSize_(I) {stencil_.resize(stencilSize_);};
        stencil(size_t I, size_t J):           I_(I), J_(J), stencilSize_(I*J) {stencil_.resize(stencilSize_);};
        stencil(size_t I, size_t J, size_t K): I_(I), J_(J), K_(K), stencilSize_(I*J*K) {stencil_.resize(stencilSize_);};

        ~stencil() {};

        void SetStencil(size_t I) {I_ = I; stencilSize_ = I; stencil_.resize(stencilSize_);}
        void SetStencil(size_t I, size_t J) {I_ = I; J_ = J; stencilSize_ = I*J; stencil_.resize(stencilSize_);}
        void SetStencil(size_t I, size_t J, size_t K) {I_ = I; J_ = J; K_ = K; stencilSize_ = I*J*K; stencil_.resize(stencilSize_);}

        void CreateStencil() {stencil_.resize(stencilSize_);}

        const int getSize() const {return stencilSize_;};

        const int getI() const {return I_;};
        const int getJ() const {return J_;};
        const int getK() const {return K_;};

        void setI(size_t I) const {I_ = I;};
        void setJ(size_t J) const {J_ = J;};
        void setK(size_t K) const {K_ = K;};

        T &operator()(size_t i){
            return stencil_[i];
        }

        T &operator()(size_t i, size_t j){
            return stencil_[i+j*I_];
        }

        T &operator()(size_t i, size_t j, size_t k){
            return stencil_[i+j*I_+k*I_*J_];
        }

        T operator()(size_t i) const{
            return stencil_[i];
        }

        T operator()(size_t i, size_t j) const{
            return stencil_[i+j*I_];
        }

        T operator()(size_t i, size_t j, size_t k) const{
            return stencil_[i+j*I_+k*I_*J_];
        }

        const vector<T>& GetStencil() const {return stencil_;}

        // Call it when T is pointer
        void clearPtr() {if (std::is_pointer<T>::value){
                             std::for_each(stencil_.begin(),stencil_.end(),delete_pointed_to<T>);}}

    private:
        size_t I_;
        size_t J_;
        size_t K_; 
        size_t stencilSize_;
        vector<T> stencil_;
};

#endif
