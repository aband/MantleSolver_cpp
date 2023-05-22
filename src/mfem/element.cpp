#include "element.h"

element::element(const MeshInfo& mi,
                 const indice& global){
    GetCorners(mi,global);
}

void element::GetCorners(const vertexSet& corners){
    for (const auto& c: corners){
        corners_.push_back(c);
    }
}

double element::lambda(const int& e,
                       const vertex& point) const{
    vertexSet edge {corners_.at(e), 
                    corners_.at((e+1)%4)};
    return distance_(edge,point);
}

double element::lambda(const int& e1,
                       const int& e2,
                       const vertex& point) const{
    

}

double element::distance_(const vertexSet& edge, 
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
void element::Test(const vertex& point){

    for (const auto& c: corners_){
        nodePrint(c);
    }cout << endl;

    for (int e = 0; e<4; e++){
        cout << "Distance is " << lambda(e,point)  << endl;
    }
}
