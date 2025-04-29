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
    corners_.clear();
    for (const auto& c: corners){
        corners_.push_back(c);
    }

    //! Create unit normal and unit tangent vectors for each edges
    unitNormals_.clear();
    unitTangents_.clear();

    for (int c=0; c<4; c++){
        vertexSet edge {corners_.at((c+3)%4), 
                        corners_.at(c)};
        unitNormals_.push_back(UnitNormal(edge, length(edge)));
        unitTangents_.push_back(unitTangent(edge, length(edge)));
    }

    unitNormals_d_.clear();

    vertexSet diag {corners_.at(0),corners_.at(2)};

    unitNormals_d_.push_back(UnitNormal(diag, length(diag)));

    diag = {corners_.at(1),corners_.at(3)};

    unitNormals_d_.push_back(UnitNormal(diag, length(diag)));

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

double basis::lambdad(const int& i,
                      const vertex& point)const{

    //! i = 0 or 1

    vertex tmp = -1*(point - corners_.at(i));

    return tmp[0]*unitNormals_d_.at(i)[0] + 
           tmp[1]*unitNormals_d_.at(i)[1]; 

}

double basis::R(const int& e1,
                const int& e2,
                const vertex& point) const {

    return (lambda(e1,point) - lambda(e2,point))/
           (lambda(e1,point) + lambda(e2,point));
}

double basis::R(const int& e,
                const vertex& point) const{
    return lambda((e+2)%4,point)/(lambda(e,point) + lambda((e+2)%4,point));
}

// Return derivative of R
vertex basis::dR(const int& e,
                 const vertex& point) const{

    vertex work(2);

    work =  (-1*lambda(e,point)*unitNormals_.at((e+2)%4) +
               unitNormals_.at(e)*lambda((e+2)%4,point)) / 
               pow(lambda(e,point)+lambda((e+2)%4,point),2);

    return work;
}

double basis::rational(const vertex& point) const{

    return lambda(2,point)*lambda(3,point)/lambda(2,corners_.at(0))/lambda(3,corners_.at(0)) - 
           lambda(3,point)*lambda(0,point)/lambda(3,corners_.at(1))/lambda(0,corners_.at(1)) + 
           lambda(0,point)*lambda(1,point)/lambda(0,corners_.at(2))/lambda(1,corners_.at(2)) -
           lambda(1,point)*lambda(2,point)/lambda(1,corners_.at(3))/lambda(2,corners_.at(3));
}

vertex basis::dRational(const vertex& point) const{

    return -1*(unitNormals_.at(2)*lambda(3,point)+
               lambda(2,point)*unitNormals_.at(3))/
           lambda(2,corners_.at(0))/lambda(3,corners_.at(0)) - 
           -1*(unitNormals_.at(3)*lambda(0,point)+
               lambda(3,point)*unitNormals_.at(0))/
           lambda(3,corners_.at(1))/lambda(0,corners_.at(1)) + 
           -1*(unitNormals_.at(0)*lambda(1,point)+
               lambda(0,point)*unitNormals_.at(1))/
           lambda(0,corners_.at(2))/lambda(1,corners_.at(2)) -
           -1*(unitNormals_.at(1)*lambda(2,point)+
               lambda(1,point)*unitNormals_.at(2))/
           lambda(1,corners_.at(3))/lambda(2,corners_.at(3));  
}

// Two lagrangian interpolation on edge nodes and vertex nodes
std::array<double, 3> basis::lagrangeE(const vertex& point, 
                                       int nEdge, int j) const{

    // First value represents value evaluated at the given point
    // Second and third values represent derivative values 
    // evaluated at the given point.
    // This lagrangian polynomial evaluates 1 on x_nEdge_j, 
    // and - on other nodes and two vertices on edge_nEdge
    std::array<double, 3> work = {1.0,0.0,0.0};

    int num_term = polynomial_degree_+1;

    std::vector<double> term_grad_coef_part(num_term,1);
    std::vector<vertex> term_grad(num_term, {0,0});

    double tmp = 0.0;

    vertex pNode = lagEdgeNode_.at(FlatIndic(polynomial_degree_-1, j, nEdge));

    // evaluation *= (point - x_{e,n,k})/(x_{e,n,j} - x_{e,n,k}) for all k != j
    for (int k=0; k<polynomial_degree_-1; k++){
        if (k == j){
            continue;
        } else {
            // Project current point onto target edge.
            double projPt = projToEdge_(nEdge, point);

            vertex zeroNode = lagEdgeNode_.at(FlatIndic(polynomial_degree_-1, k, nEdge));

            tmp = (projPt - projToEdge_(nEdge, zeroNode)) / 
                  (projToEdge_(nEdge, pNode) - projToEdge_(nEdge, zeroNode));

            work[0] *= tmp;

            for (int m=0; m<num_term; m++){
                term_grad_coef_part[m] *= (m==k) ? 1:tmp;
            }
            term_grad[k] = unitTangents_.at(nEdge) / (projToEdge_(nEdge, pNode) - 
                                                      projToEdge_(nEdge, zeroNode));
        }
    }

    // result *= (pt-x_{v,n})/(x_{e,n,j} - x_{v,n})
    for (int n=nEdge-1+4; n<nEdge+4; n++){
        vertex zeroNode = corners_.at(n%4);        
        double projPt = projToEdge_(nEdge, point);
        tmp = (projPt - projToEdge_(nEdge, zeroNode))/
              (projToEdge_(nEdge, pNode) - projToEdge_(nEdge, zeroNode));
        work[0] *= tmp;
        for (int m=0; m<num_term; m++){
            term_grad_coef_part[m] *= (m==polynomial_degree_+n-(nEdge+4)) ? 1 : tmp;
        }
        term_grad[polynomial_degree_+n-(nEdge+4)] = unitTangents_.at(nEdge) / 
                                                    (projToEdge_(nEdge,pNode) - 
                                                     projToEdge_(nEdge,zeroNode));
    }

    for (int i=0; i<num_term; i++){
        work[1] += term_grad_coef_part.at(i) * term_grad.at(i)[0];
        work[2] += term_grad_coef_part.at(i) * term_grad.at(i)[1];
    }

    return work;
}

std::array<double, 3> basis::lagrangeV(const vertex& point, 
                                       int nEdge, int i) const{

    // First value represents value evaluated at the given point
    // Second and third values represent derivative values
    // evaluated at the given point.

    assert(nEdge == i || nEdge == (i+3)%4);

    std::array<double, 3> work = {1.0,0.0,0.0};

    int num_term = polynomial_degree_;

    double tmp;

    vertex pNode = corners_.at(i);
    vertex zeroNode;

    double projPt;

    std::vector<double> term_grad_coef_part(num_term,1);
    std::vector<vertex> term_grad(num_term, {0,0});

    // result *= (pt - x_{e,n,k})/(x_{v,i} - x_{e,n,k}) for all k
    for (int k=0; k<polynomial_degree_-1; k++){
        zeroNode = lagEdgeNode_.at(FlatIndic(polynomial_degree_-1, k, nEdge));
        projPt = projToEdge_(nEdge, point);

        tmp = (projPt - projToEdge_(nEdge, zeroNode)) / 
              (projToEdge_(nEdge, pNode) - projToEdge_(nEdge, zeroNode));
        work[0] *= tmp;

        for (int n=0; n<num_term; n++){
            term_grad_coef_part[n] *= (n==k) ? 1 : tmp;
        }
        term_grad[k] = unitTangents_.at(nEdge) / (projToEdge_(nEdge, pNode) - 
                                                  projToEdge_(nEdge, zeroNode));
    }

    if (nEdge == i){
        zeroNode = corners_.at((i+3)%4);
    } else {
        zeroNode = corners_.at((i+1)%4);
    }

    projPt = projToEdge_(nEdge, point);

    tmp = (projPt - projToEdge_(nEdge, zeroNode))/
          (projToEdge_(nEdge, pNode) - projToEdge_(nEdge, zeroNode));

    work[0] *= tmp; 

    for (int n=0; n<num_term; n++){
        term_grad_coef_part[n] *= (n==polynomial_degree_-1) ? 1:tmp;
    }

    term_grad[polynomial_degree_-1] = unitTangents_.at(nEdge) / 
                                      (projToEdge_(nEdge,pNode) - projToEdge_(nEdge, zeroNode));

    for (int i=0; i<num_term; i++){
        work[1] += term_grad_coef_part.at(i) * term_grad.at(i)[0];
        work[2] += term_grad_coef_part.at(i) * term_grad.at(i)[1];
    }

    return work;
}

// Supplemental functions
double basis::PhiSupp(const int& i,
                      const vertex& point,
                      const int& r) const{

    assert(r>2);

    switch(i){
        case 0:
            return lambda(2-1,point)*lambda(4-1,point)*pow(lambda(2-1,4-1,point),r-2)*R(1,3,point);
        case 1:
            return lambda(1-1,point)*lambda(3-1,point)*pow(lambda(1-1,3-1,point),r-2)*R(0,2,point);
        default:
            cout << "Supplement function not defined. " << endl;
            return -1;
    }

}

// =======================================================================
double basis::projToEdge_(int Edge, const vertex& point) const {

    // The derivative of this projection is just tau 

    vertex tmp = point - corners_.at((Edge+3)%4); 

    return std::inner_product(std::begin(tmp),
                              std::end(tmp),
                              std::begin(unitTangents_.at(Edge)),
                              0.0);   
}


double basis::distance_(const vertexSet& edge, 
                        const vertex& point) const{

    //vertex tmp = point - (edge.at(0) + edge.at(1))/2;
    vertex tmp = point - edge.at(1);

    vertex unitNormal = UnitNormal(edge, length(edge));

    //return -1*std::inner_product(std::begin(tmp),
    //                             std::end(tmp),
    //                             std::begin(unitNormal),
    //                             0.0);

    return -1*(unitNormal[0]*tmp[0]+unitNormal[1]*tmp[1]);
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

    cout << "Test projection of points" << endl;

    vertex tmp {0,0};

    double projP = projToEdge_(0, tmp); 
    cout << projP << endl;
    projP = projToEdge_(1, tmp); 
    cout << projP << endl;
    projP = projToEdge_(2, tmp); 
    cout << projP << endl;
    projP = projToEdge_(3, tmp); 
    cout << projP << endl;
}
