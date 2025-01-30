#include "tensor.h"

// Include some tensor arrithmatic

int Tensor_add(const Tensor<double>& t1,
               const Tensor<double>& t2,
               Tensor<double>& t3){

    assert(t1.getSize() == t2.getSize());
    assert(t1.getSize() == t3.getSize());

    for (int i=0; i<t3.getSize(); i++){
        t3(i) = t1(i) + t2(i);
    }

    return 1;
}

int Tensor_multi(const Tensor<double>& t1, 
                 const Tensor<double>& t2, 
                 Tensor<double>& t3){

    assert(t1.getSize() == t2.getSize());
    assert(t1.getSize() == t3.getSize());

    for (int i=0; i<t3.getSize(); i++){
        t3(i) = t1(i) * t2(i);
    }

    return 1;
}

int Tensor_zero(Tensor<double>& t){

    for (int i=0;i<t.getSize(); i++){
        t(i) = 0.0;
    }

    return 1;
}

int Tensor_multi_add(const Tensor<double>& t1, 
                     const Tensor<double>& t2,
                     double scale,
                     Tensor<double>& t3){

    assert(t1.getSize() == t2.getSize());
    assert(t1.getSize() == t3.getSize());

    for (int i=0; i<t3.getSize(); i++){
        t3(i) += t1(i) * t2(i) * scale;
    }

    return 1;
}

int Tensor_scale(double scale, Tensor<double>& t){

    for (int i=0; i<t.getSize(); i++){
        t(i) *= scale;
    }

    return 1;
}
