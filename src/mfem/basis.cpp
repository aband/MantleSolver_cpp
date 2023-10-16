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
    vertexSet edge {corners_.at((e+3)%4), 
                    corners_.at(e)};
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
    return 0.5*(1-R(e,(e+2)%4,point));
}

// Supplemental functions
double basis::PhiSupp(const int& i,
                      const vertex& point,
                      const int& r) const{

    assert(r>2);

    switch(i){
        case 0:
            return lambda(2-1,point)*lambda(4-1,point)*pow(lambda(2-1,4-1),r-2,point)*R(1,3,point);
        case 1:
            return lambda(1-1,point)*lambda(3-1,point)*pow(lambda(1-1,3-1),r-2,point)*R(0,2,point);
        default:
            cout << "Supplement function not defined. " << endl;
            return -1;
    }

}

// Edge nodal basis functions
double basis::phi_e() const{

    double work = 0.0;



    return work;
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
        cout << "Distance is " << lambda(e,point)  << endl;
    }

    // Test R function
    int seed = 10;
    std::array<int,2> list {0, 2};

    for (int Case = 0; Case<4; Case ++){
        cout << "Current edge is: " << Case << endl;
        vertexSet edge {corners_.at((Case+3)%4), corners_.at(Case)};
        vertex unittangent = unitTangent(edge, length(edge));
        double dl = length(edge) / seed;
        //cout << "On the edge " << e+1 << " the values are distributed as: " ;
        for (int j=0; j<seed; j++){
            cout << R(Case,corners_.at(Case) + j*unittangent*dl) << " "; 
        } cout << endl;

        // ======================================================================
        cout << "Opposite edge is: " << (Case+2)%4 << endl;
        edge = {corners_.at((Case+3+2)%4), corners_.at((Case+2)%4)};
        unittangent = unitTangent(edge, length(edge));
        dl = length(edge) / seed;
        //cout << "On the edge " << e+1 << " the values are distributed as: " ;
        for (int j=0; j<seed; j++){
            cout << R(Case,corners_.at((Case+2)%4) + j*unittangent*dl) << " "; 
        } cout << endl;
    }
}
