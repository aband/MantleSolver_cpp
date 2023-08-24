#include "interp.h"

interpolation::Lagrange::Lagrange(const tensor<int>& order){
    for (int it = 0; it<order.size(); it++){
        order_[it] = order[it];
    }
}

void interpolation::Lagrange::init(const tensorSet<double>& points){
    assert(interpoPoints_.size() == 0);

    interpoPoints_.resize(points.size());

    for (auto it : points){
        interpoPoints_.push_back(it);
    }
}

void interpolation::Lagrange::clear(){
    interpoPoints_.clear();
}

vertexSet interpolation::Lagrange::chebyshevPoints_(const vertex& left, 
                                                    const vertex& right){
    // Create chebyshev points on a given edge
    // The edge is defined by left and right vertices

    vertexSet points;


    return points;
}

double interpolation::Lagrange::eval(const tensor<double>& point) const {

    double work = 0.0;

    // Barycentric lagrange interpolatin

    return work;
}

double interpolation::Lagrange::evalDeriv(const tensor<double>& point) const {

    double work = 0.0;

    return work;
}
