#include "basis.h"

//=================================
//       e4 
//   v4 ----- v3
//   |        |
// e1|        | e3
//   |        |
//   v1 ----- v2
//       e2
//=================================

// ===== Public members =====
void basis::GetCorners(const vertexSet& corners){
    for (const auto& c: corners){
        corners_.push_back(c);
    }
}

double basis::lambda(const int& e,
                     const vertex& point) const{
    // e = {1,2,3,4}
    assert(e>0 & e<5);
    vertexSet edge {corners_.at((e+2)%4), 
                    corners_.at((e+3)%4)};
    return distance_(edge,point);
}

double basis::lambda(const int& e1,
                     const int& e2,
                     const vertex& point) const{

    vertexSet line {corners_.at(e1),
                    corners_.at(e2)};

    return (lambda(e1,point) - lambda(e2,point))/length(line);
}

double basis::R(const int& e1,
                const int& e2,
                const vertex& point) const {

    return (lambda(e1,point) - lambda(e2,point))/
           (lambda(e1,point) + lambda(e2,point));
}

double basis::R(const int& e,
                const vertex& point) const{
    return 0.5*(1-R(e,(e+3)%4,point));
}

// Supplemental functions
double basis::PhiSupp(const int& i,
                      const vertex& point) const{

    return lambda(i,point)*lambda(i+2,point);
}


double basis::phi() const{

    double work;

    return work;
}

// ===== Private members =====
double basis::distance_(const vertexSet& edge, 
                        const vertex& point) const{

    //vertex tmp = point - (edge.at(0) + edge.at(1))/2;
    vertex tmp = point - edge.at(0);

    vertex unitNormal = UnitNormal(edge, length(edge));

    return -1*std::inner_product(std::begin(tmp),
                                 std::end(tmp),
                                 std::begin(unitNormal),
                                 0.0);
}

// ===== Test =====
void basis::Test(const vertex& point){

    for (const auto& c: corners_){
        nodePrint(c);
    }cout << endl;

    for (int e = 0; e<4; e++){
        cout << "Distance is " << lambda(e+1,point)  << endl;
    }
}
