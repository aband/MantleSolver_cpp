#include "shape.h"

template <typename T> std::vector<vertex> 
shape<T>::EvaluateAll(const vertex& point) const{

    return ptr_->EvaluateAll(*basisPtr_, point);
}

template <typename T> 
vertex shape<T>::Evaluate(const vertex& point, 
                          const int& localdof) const{

    return ptr_->Evaluate(*basisPtr_, point, localdof);
}

template <typename T> 
std::vector<vertex> shape<T>::EvaluateGradAll(const vertex& point) const{
		  
    return ptr_->EvaluateGradAll(*basisPtr_, point);
}

template <typename T> 
std::vector<int> shape<T>::LocalToGlobal(const MeshInfo& mi,
                                         const indice& global) const{

    return ptr_->LocalGlobalMap(mi, global);
}

template <typename T>
std::vector<int> shape<T>::GlobalToLocalMapBndry(const MeshInfo& mi,
                                                 const int& gdof) const{

    return ptr_->GlobalToLocalMapBndry(mi, gdof);
}
