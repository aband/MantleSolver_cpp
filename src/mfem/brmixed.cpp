#include "brmixed.h"

// Local Degree of Freedom
// x, y, supplement

//
//  3_____2    7_____6       _11__
//  |     |    |     |      |     | 
//  |     |    |     |     8|     | 10
//  0-----1    4-----5       -----
//                             9

// Return a global index of the given 
std::array<int,12> BRMixed::LocalToGlobal(const MeshInfo& mi,
                                          const indice& globalElement) const{

    std::array<int, 12> globalDOF;

    globalDOF[0] = FlatIndic(mi.MPIglobalVertexSize[0],globalElement); 
    globalDOF[1] = globalDOF[0] + 1;

    indice add {1,1};

    indice tmp = globalElement + add;

    globalDOF[2] = FlatIndic(mi.MPIglobalVertexSize[0],tmp);
    globalDOF[3] = globalDOF[2] - 1;

    int totalNodal = mi.MPIglobalVertexSize[0] * mi.MPIglobalVertexSize[1];

    globalDOF[4] = totalNodal + globalDOF[0];
    globalDOF[5] = totalNodal + globalDOF[1];
    globalDOF[6] = totalNodal + globalDOF[2];
    globalDOF[7] = totalNodal + globalDOF[3];

    globalDOF[8] = totalNodal*2 + mi.MPIglobalHoriEdgeSize + 
                   FlatIndic(mi.MPIglobalVertexSize[0],globalElement);

    globalDOF[9] = totalNodal*2 + FlatIndic(mi,globalElement);

    globalDOF[10] = globalDOF[8] + 1;

    add = {0,1};

    tmp = globalElement + add;

    globalDOF[11] = totalNodal*2 + FlatIndic(mi,tmp);

    return globalDOF;
}

void BRMixed::ComputeTotalDOF(const MeshInfo& mi){
    totalDOF_ = 2*mi.MPIglobalVertexSize[0] * mi.MPIglobalVertexSize[1] + 
                mi.MPIglobalHoriEdgeSize + mi.MPIglobalVertEdgeSize;
}

double BRMixed::phiv(const basis& basis_,
                     const int& nnodal,
                     const vertex& point) const {

    int idiag = (nnodal+1)%2;

    double work = 0.0;

    work = basis_.lambdad(idiag,point) - 
       0.5*basis_.lambdad(idiag,basis_.corners_.at((nnodal+2)%4))*
       (1+pow(-1,idiag+1) * R_(basis_,point)); 

    work /= basis_.lambdad(idiag,basis_.corners_.at(nnodal)) - 
            basis_.lambdad(idiag,basis_.corners_.at((nnodal+2)%4));

    return work;

}

double BRMixed::phie(const basis& basis_,
                     const int& nEdge,
                     const vertex& point) const {

    vertex mid = (basis_.corners_.at((nEdge+3)%4) + basis_.corners_.at(nEdge))/2.0;

    return phie_(basis_,nEdge,point)/phie_(basis_,nEdge,mid);
}

vertex BRMixed::dPhie(const basis& basis_,
                      const int& nEdge,
                      const vertex& point) const{

    vertex mid = (basis_.corners_.at((nEdge+3)%4) + basis_.corners_.at(nEdge))/2.0;

    return dphie_(basis_,nEdge,point)/phie_(basis_,nEdge,mid);

}

// Gradient of a given vector
// grad[a] = [a_x, a_y]
//     [b]   [b_x, b_y]

vertex BRMixed::dPhiv(const basis& basis_,
                      const int& nnodal,
                      const vertex& point) const{

    // Compute derivative of the given basis function
    // [f_x,f_y] first value is x derivative and the second value is y derivative

    int idiag = (nnodal+1)%2;

    vertex work(2);

    work = -1*basis_.unitNormals_d_.at(idiag) - 
          0.5*basis_.lambdad(idiag,basis_.corners_.at((nnodal+2)%4))*
             (pow(-1,idiag+1) * dR_(basis_,point)); 

    work /= basis_.lambdad(idiag,basis_.corners_.at(nnodal)) - 
            basis_.lambdad(idiag,basis_.corners_.at((nnodal+2)%4));

    return work;
}

// Calculate gradient of corresponding basis functions
std::array<std::array<double,4>, 12> BRMixed::ComputeGradBRmixed(const basis& basis_,
                                                                 const vertex& point)const {

    std::array<std::array<double,4>, 12> work;

    for (unsigned int j=0; j<4; j++){
        vertex dphiv = dPhiv(basis_,j,point); 

        work[j][0] = dphiv[0];
        work[j][1] = dphiv[1];
        work[j][2] = 0.0;
        work[j][3] = 0.0;

        work[j+4][0] = 0.0;
        work[j+4][1] = 0.0;
        work[j+4][2] = dphiv[0];
        work[j+4][3] = dphiv[1];

        vertex unitN = basis_.unitnormal(j);
        vertex dphie = dPhie(basis_,j,point); 

        work[j+8][0] = dphie[0] * unitN[0];
        work[j+8][1] = dphie[1] * unitN[0];
        work[j+8][2] = dphie[0] * unitN[1];
        work[j+8][3] = dphie[1] * unitN[1];

        //cout << work[j+8][0] << "  " << work[j+8][1] << "  " << work[j+8][2] << "  " << work[j+9][3] << endl; 

    }

    // Fix global unit normal edge directions
    for (unsigned int k=0; k<4; k++){
        work[8][k] *= -1;
        work[9][k] *= -1;
    }

    return work;
}

std::array<vertex, 12> BRMixed::ComputeBRmixed(const basis& basis_,
                                               const vertex& point) const{

    std::array<vertex, 12> work;

    for (unsigned int k=0; k<4; k++){
        work[k]   = {phiv(basis_,k,point), 0.0};
        work[k+4] = {0.0, phiv(basis_,k,point)};
        work[k+8] = basis_.unitnormal(k) * phie(basis_,k,point); 
    }

    work[8] *= -1;
    work[9] *= -1;

    return work;
}

vertex BRMixed::ComputeBRmixed(const basis& basis_,
                               const vertex& point,
                               const int& local) const{

    vertex work;

    std::array<double, 4> coef {-1.0,-1.0,0.0,0.0};

    if (local < 4){
        work = {phiv(basis_,local,point), 0.0};
    } else if (local < 8){
        work = {0.0, phiv(basis_,local-4,point)};
    } else {
        work = basis_.unitnormal(local-8) * 
               phie(basis_,local-8,point) * coef[local-8];
    }

    return work;
}

// =============================================================================

std::vector<vertex> BRMixed::EvaluateAll(const basis& basis_,
                                         const vertex& point) const{
    std::array<vertex, 12> tmp = ComputeBRmixed(basis_, point);
    std::vector<vertex> work (tmp.begin(), tmp.end()); 

    return work;
}

vertex BRMixed::Evaluate(const basis& basis_,
                         const vertex& point, 
                         const int& localdof) const{

    return ComputeBRmixed(basis_,point,localdof);
}

std::vector<std::array<double,4>> BRMixed::EvaluateGradAll(const basis& basis_,
                                                           const vertex& point) const{

    std::array<std::array<double, 4>, 12> tmp = ComputeGradBRmixed(basis_, point); 

    std::vector<std::array<double,4>> work (tmp.begin(), tmp.end());

    return work;
}

std::vector<int> BRMixed::LocalGlobalMap(const MeshInfo& mi,
                                         const indice& global) const{
  
    std::array<int, 12> tmp = LocalToGlobal(mi, global);
    std::vector<int> work (tmp.begin(), tmp.end());
    return work;
}

std::vector<int> BRMixed::GlobalToLocalMapBndry(const MeshInfo& mi,
                                                const int& globaldof) const{

    //! Get Cell index and corresponding local dof index
    //! Will only return cell index and local dof for those dofs right on boundary
    //! Attention !!!!
    //! Interior dofs will return {-1} not cell index and local dof

    std::vector<int> work;

    int totalNodal = mi.MPIglobalVertexSize[0] * mi.MPIglobalVertexSize[1];

    int moddof;

    indice bend;

    if (globaldof > totalNodal*2-1){
        // It is a bubble function dof
        moddof = globaldof - totalNodal*2;

        if (moddof > mi.MPIglobalHoriEdgeSize-1){

            moddof -= mi.MPIglobalHoriEdgeSize;
            bend = Bend(mi.MPIglobalVertexSize[0],moddof); 

            if (bend[0] == 0){
                // On boundary dof.
                // left 
                work.push_back(FlatIndic(mi,bend));
                work.push_back(8);

            } else if (bend[0] == mi.MPIglobalVertexSize[0]-1){
                // On boundary dof
                // right
                work.push_back(FlatIndic(mi,bend[0]-1,bend[1]));
                work.push_back(10);

            } else {
                // Interior dof
                work.push_back(-1);
            }

        } else {
            bend = Bend(mi.MPIglobalCellSize[0], moddof);
            if (bend[1] == 0 ){
                // On boundary dof.
                // bottom 
                work.push_back(FlatIndic(mi,bend));
                work.push_back(9);

            } else if (bend[1] == mi.MPIglobalVertexSize[1]-1){
                // On boundary dof.
                // top
                work.push_back(FlatIndic(mi,bend[0],bend[1]-1));
                work.push_back(11);
               
            } else {
                work.push_back(-1);
            }
        }

    } else {

        int shift = 0;
        if (globaldof > totalNodal-1){
            // It is a y direction dof
            moddof = globaldof - totalNodal;
				shift  = 4;
        } else {
            // It is a x direction dof
            moddof = globaldof;
            shift  = 0;
        }
        // x and y dofs are treated similarly

        bend = Bend(mi.MPIglobalVertexSize[0], moddof);

        if (bend[0] == 0 && 
            bend[1] != 0){
            // top left
            work.push_back(FlatIndic(mi,bend[0],bend[1]-1));
            work.push_back(3+shift); 

        } else if(bend[1] == 0 &&
                  bend[0] != mi.MPIglobalVertexSize[1]-1){
            // bottom left
            work.push_back(FlatIndic(mi,bend[0],bend[1]));
            work.push_back(0+shift); 

        } else if(bend[0] == mi.MPIglobalVertexSize[0]-1 && 
                  bend[1] != mi.MPIglobalVertexSize[1]-1){
            // bottom right
            work.push_back(FlatIndic(mi,bend[0]-1,bend[1]));
            work.push_back(1+shift); 

        } else if(bend[1] == mi.MPIglobalVertexSize[1]-1 &&
                  bend[0] != 0){
            // top right
            work.push_back(FlatIndic(mi,bend[0]-1,bend[1]-1));
            work.push_back(2+shift);

        } else {
            work.push_back(-1);
        }

        // manually fixing right bottom corner
        if (bend[0] == mi.MPIglobalVertexSize[0] -1 &&
            bend[1] == 0){
            work[0] = mi.MPIglobalCellSize[0]-1;
            work[1] = 1+shift;
        } 

        if (bend[0] == 0 &&
            bend[1] == 1){
            work[0] = mi.MPIglobalCellSize[0];
            work[1] = 0+shift;
        }

    }

    return work;
}

bool BRMixed::onBndry(const MeshInfo& mi,
                      const int& globaldof) const{

    bool result = false;

    int totalNodal = mi.MPIglobalVertexSize[0] * mi.MPIglobalVertexSize[1];

    int moddof;

    indice bend;

    if (globaldof > totalNodal*2-1){
        // It is a bubble function dof
        moddof = globaldof - totalNodal*2;

        // Count horizontal edges first
        if (moddof > mi.MPIglobalHoriEdgeSize-1){
            moddof -= mi.MPIglobalHoriEdgeSize;
            bend = Bend(mi.MPIglobalVertexSize[0],moddof); 

            if (bend[0] == 0 || bend[0] == mi.MPIglobalVertexSize[0]-1){
                result = true;
            }

        } else {
            bend = Bend(mi.MPIglobalCellSize[0], moddof);
            if (bend[1] == 0 || bend[1] == mi.MPIglobalVertexSize[1]-1){
                result = true;
            }
        }

    }else if (globaldof > totalNodal-1){
        moddof = globaldof - totalNodal;
        bend = Bend(mi.MPIglobalVertexSize[0], moddof);

        if (bend[0] == 0 || bend[1] == 0 || 
            bend[0] == mi.MPIglobalVertexSize[0]-1 ||
            bend[1] == mi.MPIglobalVertexSize[1]-1){
            result = true;
        }

    }else {

        moddof = globaldof;
        bend = Bend(mi.MPIglobalVertexSize[0], moddof);

        if (bend[0] == 0 || bend[1] == 0 || 
            bend[0] == mi.MPIglobalVertexSize[0]-1 ||
            bend[1] == mi.MPIglobalVertexSize[1]-1){
            result = true;
        }

    }

    return result;
}

// ==================================================================

double BRMixed::phie_(const basis& basis_,
                      const int& nEdge,
                      const vertex& point) const{

    return basis_.lambda((nEdge+1)%4,point)*basis_.lambda((nEdge+3)%4,point)*
           basis_.R(nEdge,point);
}

double BRMixed::R_(const basis& basis_,
                   const vertex& point) const {

    double work = basis_.rational(point);

    for (int nedge = 0; nedge<4; nedge ++){
        vertex mid = (basis_.corners_.at((nedge+3)%4) + 
                      basis_.corners_.at(nedge))/2.0;

        work -= basis_.rational(mid)*phie(basis_,nedge,point);
    }

    return work;
}

vertex BRMixed::dphie_(const basis& basis_,
                       const int& nEdge,
                       const vertex& point) const {
    return  -1*basis_.unitNormals_.at((nEdge+1)%4)*basis_.lambda((nEdge+3)%4,point) * basis_.R(nEdge,point) -
               basis_.unitNormals_.at((nEdge+3)%4)*basis_.lambda((nEdge+1)%4,point) * basis_.R(nEdge,point) +
               basis_.lambda((nEdge+1)%4,point) * basis_.lambda((nEdge+3)%4,point) * basis_.dR(nEdge,point);

/*
    return (-1*(basis_.lambda(nEdge,point)+basis_.lambda((nEdge+2)%4,point)) *
           (basis_.lambda((nEdge+1)%4,point) * basis_.lambda((nEdge+3)%4,point) * basis_.unitNormals_.at((nEdge+2)%4) + 
            basis_.lambda((nEdge+1)%4,point) * basis_.lambda((nEdge+2)%4,point) * basis_.unitNormals_.at((nEdge+3)%4) + 
            basis_.lambda((nEdge+2)%4,point) * basis_.lambda((nEdge+3)%4,point) * basis_.unitNormals_.at((nEdge+1)%4)) + 
            basis_.lambda((nEdge+2)%4,point) * basis_.lambda((nEdge+3)%4,point) * basis_.lambda((nEdge+4)%4,point) * (basis_.unitNormals_.at((nEdge+0)%4) + basis_.unitNormals_.at((nEdge+2)%4)))/ 
            pow(basis_.lambda(nEdge,point) + basis_.lambda((nEdge+2)%4,point),2);
*/

}

vertex BRMixed::dR_(const basis& basis_,
                    const vertex& point) const{
    vertex work = basis_.dRational(point);

    for (int nedge = 0; nedge<4; nedge++){
        vertex mid = (basis_.corners_.at((nedge+3)%4) + 
                      basis_.corners_.at(nedge))/2.0;

        work -= basis_.rational(mid)*dPhie(basis_,nedge,point);

    }

    return work;
}

void BRMixed::Test(const basis& basis_){

    int seed = 10;

    for (int Case = 0; Case<4; Case++){
        cout << "Current edge is: " << Case << endl;

        vertexSet edge {basis_.corners_.at((Case+3)%4), basis_.corners_.at(Case)};
        vertex unittangent = unitTangent(edge, length(edge));
        double len = length(edge);
        double dl = len / seed;

        for (int j=0; j<seed+1; j++){
            vertex point =  basis_.corners_.at((Case+3)%4) + j*unittangent*dl;

            vertex mid = 0.5*(basis_.corners_.at((Case+3)%4) + 
                              basis_.corners_.at(Case));

            double val = phiv(basis_,Case, point);

            cout << val << " " ;

        } cout << endl;

    }
}
