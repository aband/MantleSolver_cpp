#ifndef SHAPE_H_
#define SHAPE_H_

#include "basis.h"
#include "Hdivmixed.h"
#include "brmixed.h"

template <typename T> 
class shape{

    public :
       shape(basis * basisptr, T * ptr) {basisPtr_ = basisptr; ptr_ = ptr;};
       ~shape() {};

       // ! Standarized evaluation functions
       std::vector<vertex> EvaluateAll(const vertex& point) const;

       vertex Evaluate(const vertex& point, 
                       const int& localdof) const;

       std::vector<vertex> EvaluateGradAll(const vertex& point) const;

       std::vector<int> LocalToGlobal(const MeshInfo& mi,
                                      const indice& global) const;

       int getDOF() const {return ptr_->getDOF();}; 

       std::vector<int> GlobalToLocalMapBndry(const MeshInfo& mi,
                                              const int& gdof)const;

       bool onBndry(const MeshInfo& mi, 
            const int& globaldof) const {return ptr_->onBndry(mi, globaldof);};

       double Pressure() const {return ptr_.Pressure();};

       vertexSet corners() const {return basisPtr_->corners();};

       std::string name() const{return ptr_->name;}

    private:
       // Pointer to the shape function in use    
       T * ptr_;

       basis * basisPtr_;
};

#include "shape.tpp"

#endif
