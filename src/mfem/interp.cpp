#include "interp.h"

Lagrange::interpolation::interpolation(const tensor<int>& order){
    for (int it = 0; it<order.size(); it++){
        order_[it] = order[it];
    }
}

void Lagrange::interpolation::init(const tensorSet<double>& points){
    assert(interpoPoints_.size() == 0);

    interpoPoints_.resize(points.size());

    for (auto it : points){
        interpoPoints_.push_back(it);
    }
}

void Lagrange::interpolation::clear(){
    interpoPoints_.clear();
}

void Lagrange::interpolation::chebyshevPoints_(){


}

double Lagrange::interpolation::eval(const tensor<double>& point) const {

    double work = 0.0;

    // Barycentric lagrange interpolatin

    return work;
}

double Lagrange::interpolation::evalDeriv(const tensor<double>& point) const {

    double work = 0.0;

    return work;
}
