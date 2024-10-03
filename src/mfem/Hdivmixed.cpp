#include "Hdivmixed.h"

std::array<int, 8> Hdivmixed::LocalToGlobal(const MeshInfo& mi,
                                            const indice& globalElement)const{

    std::array<int, 8> work;

    work[0] = mi.MPIglobalHoriEdgeSize + 
              FlatIndic(mi.MPIglobalVertexSize[0],globalElement);

    work[1] = FlatIndic(mi,globalElement);

    work[2] = work[0] + 1;

    indice add {0,1};

    indice tmp = globalElement + add;

    work[3] = FlatIndic(mi.MPIglobalCellSize[0],tmp);

    int allEdge = mi.MPIglobalHoriEdgeSize + mi.MPIglobalVertEdgeSize;

    work[4] = allEdge + work[0];

    work[5] = allEdge + work[1];

    work[6] = allEdge + work[2];

    work[7] = allEdge + work[3];

    return work;
}

void Hdivmixed::ComputeTotalDOF(const MeshInfo& mi){

    totalDOF_ = 2*(mi.MPIglobalHoriEdgeSize + mi.MPIglobalVertEdgeSize);
}

vertex Hdivmixed::phic(const basis& basis_,
                       const int& nEdge,
                       const vertex& point) const {

    vertexSet tmpEdge {basis_.corners_.at((nEdge+3)%4), basis_.corners_.at(nEdge)};

    double len1 = length(tmpEdge);

    tmpEdge = {basis_.corners_.at((nEdge+4)%4), basis_.corners_.at((nEdge+1)%4)};

    double len2 = length(tmpEdge);

    vertex work;

    vertex tmp = basis_.corners_.at(nEdge) - basis_.corners_.at((nEdge+2)%4);

    double tmpval1 = tmp[0]*basis_.unitNormals_.at(nEdge)[0] + 
                     tmp[1]*basis_.unitNormals_.at(nEdge)[1];

    double tmpval2 = tmp[0]*basis_.unitNormals_.at((nEdge+1)%4)[0] + 
                     tmp[1]*basis_.unitNormals_.at((nEdge+1)%4)[1];


    work = tmpval2 * len2 * curlphiv_star_(basis_,nEdge,point) + 
           phiv_star_star_(basis_,nEdge,point);

    work = work/(tmpval2*len2/len1 + tmpval1);

    return work;
}

vertex Hdivmixed::phil_(const basis& basis_,
                        const int& nEdge, 
                        const vertex& point) const {

    return (curlLambda_(basis_,(nEdge+1)%4)*basis_.lambda((nEdge+3)%4,point) + 
            curlLambda_(basis_,(nEdge+3)%4)*basis_.lambda((nEdge+1)%4,point))*
           basis_.R(nEdge,point) + 
           basis_.lambda((nEdge+1)%4,point)*basis_.lambda((nEdge+3)%4,point)*
           curlR_(basis_,nEdge,point);
}

vertex Hdivmixed::phil(const basis& basis_,
                       const int& nEdge,
                       const vertex& point) const {

    // Returns a normalized result

    vertexSet edge {basis_.corners_.at((nEdge+3)%4), basis_.corners_.at(nEdge)};

    vertex mid = (edge[0]+edge[1])/2.0;

    double len = length(edge); 

    return phil_(basis_, nEdge, point) / phie_(basis_, nEdge, mid) * len /4;

}

double Hdivmixed::divphic(const basis& basis_,
                          const int& nEdge,
                          const vertex& point) const{

    vertexSet tmpEdge {basis_.corners_.at((nEdge+3)%4), basis_.corners_.at(nEdge)};

    double len1 = length(tmpEdge);

    tmpEdge = {basis_.corners_.at((nEdge+4)%4), basis_.corners_.at((nEdge+1)%4)};

    double len2 = length(tmpEdge);

    double work = 0.0;

    vertex tmp = basis_.corners_.at(nEdge) - basis_.corners_.at((nEdge+2)%4);

    double tmpval1 = tmp[0]*basis_.unitNormals_.at(nEdge)[0] + 
                     tmp[1]*basis_.unitNormals_.at(nEdge)[1];

    double tmpval2 = tmp[0]*basis_.unitNormals_.at((nEdge+1)%4)[0] + 
                     tmp[1]*basis_.unitNormals_.at((nEdge+1)%4)[1];


    work = 2;

    work = work/(tmpval2*len2/len1 + tmpval1);

    return work;
}

std::array<vertex,8> Hdivmixed::ComputeHdivmixed(const basis& basis_,
                                                 const vertex& point) const{

    std::array<vertex, 8> work;

    for (unsigned int k=0; k<4; k++){
        vertex vlinear = phil(basis_, k, point);
        vertex vconst = phic(basis_, k, point);

        work[k] = vlinear;

        work[k+4] = vconst;

        /*
        work[k] = (vlinear+vconst)/2;

        work[k+4] = (vlinear-vconst)/2;
        */

    }

    // Control with global unit normal direction

    for (unsigned int k=0; k<2; k++){
        work[k] *= 1;

        work[k+4] *= -1;
    }


    return work;
}

std::array<vertex,2> Hdivmixed::ComputeHdivmixed(const basis& basis_,
                                                 const vertex& point,
                                                 const int& edge) const{

    std::array<vertex, 2> work;


    if(edge < 2){
        // Corrected with unit normal direction
        work[0] = 1*phil(basis_,edge,point);
        work[1] = -1*phic(basis_,(edge+4)%4,point);
    } else {
        work[0] = phil(basis_,edge,point);
        work[1] = phic(basis_,(edge+4)%4,point);
    }

    return work;

}

// ========================================================================================

std::vector<vertex> Hdivmixed::EvaluateAll(const basis& basis_,
                                           const vertex& point)const{

    std::array<vertex,8> tmp = ComputeHdivmixed(basis_, point);
    std::vector<vertex>  work (tmp.begin(), tmp.end());

    return work;
}

vertex Hdivmixed::Evaluate(const basis& basis_,
                           const vertex& point,
                           const int& localdof)const {

    std::array<vertex,8> tmp = ComputeHdivmixed(basis_, point);

    return tmp.at(localdof); 
} 

std::vector<std::array<double,4>> Hdivmixed::EvaluateGradAll(const basis& basis_,
                                                             const vertex& point)const{
    std::vector<std::array<double,4>> work;

    cout << "Evaluate gradient of Hdiv space not defined. Incorrectly called \n" << endl;

    return work;
}

std::vector<int> Hdivmixed::LocalGlobalMap(const MeshInfo& mi,
                                           const indice& global) const{

    std::array<int, 8> tmp = LocalToGlobal(mi, global);
    std::vector<int> work (tmp.begin(), tmp.end());
    return work;
}

std::vector<int> Hdivmixed::GlobalToLocalMapBndry(const MeshInfo& mi,
                                                  const int& globaldof)const{

    std::vector<int> work;

    int allEdge = mi.MPIglobalHoriEdgeSize + mi.MPIglobalVertEdgeSize;

    int moddof = 0;

    int shift  = 0;

    if (globaldof > allEdge-1) {
        // This is the second dof on edge
        moddof = globaldof - allEdge;
        shift  = 4;
    } else {
        moddof = globaldof;
        shift  = 0;
    }

    // Horizontal edges are counted first
    if (moddof > mi.MPIglobalHoriEdgeSize - 1){
        // It is a dof on vetical edge
        moddof -= mi.MPIglobalHoriEdgeSize;
        indice bend = Bend(mi.MPIglobalVertexSize[0], moddof);

        if (bend[0] == 0){

            work.push_back(FlatIndic(mi, bend[0], bend[1]));
            work.push_back(0+shift);

        } else if (bend[0] == mi.MPIglobalVertexSize[0] -1){
            work.push_back(FlatIndic(mi, bend[0]-1, bend[1]));
            work.push_back(2+shift);

        } else {
            // Interior
            work.push_back(-1);
        }

    } else {
        // It is a dof on horizontal edge
        indice bend = Bend(mi.MPIglobalCellSize[0], moddof);
        if (bend[1] == 0){

            work.push_back(FlatIndic(mi, bend[0],bend[1]));
            work.push_back(1+shift);

        } else if (bend[1] == mi.MPIglobalVertexSize[1]-1){

            work.push_back(FlatIndic(mi, bend[0], bend[1]-1));
            work.push_back(3+shift);

        } else {
            // Interior
            work.push_back(-1);
        }
    }

    return work;
}

bool Hdivmixed::onBndry(const MeshInfo& mi,
                        const int& globaldof) const{

    // return a bool variable determining whether this global dof is 
    // on boundary or not.

    bool result = false;

    int allEdge = mi.MPIglobalHoriEdgeSize + mi.MPIglobalVertEdgeSize;

    int moddof = 0;

    if (globaldof > allEdge-1) {
        // This is the second dof on edge
        moddof = globaldof - allEdge;
    } else {
        moddof = globaldof;
    }

    // Horizontal edges are counted first
    if (moddof > mi.MPIglobalHoriEdgeSize - 1){
        // It is a dof on vetical edge
        moddof -= mi.MPIglobalHoriEdgeSize;
        indice bend = Bend(mi.MPIglobalVertexSize[0], moddof);

        if (bend[0] == 0 || bend[0] == mi.MPIglobalVertexSize[0]-1){
            result = true;
        }

    } else {
        // It is a dof on horizontal edge
        indice bend = Bend(mi.MPIglobalCellSize[0], moddof);
        if (bend[1] == 0 || bend[1] == mi.MPIglobalVertexSize[1]-1){
            result = true;
        }
    }

    return result;
}

// ===============================================================

vertex Hdivmixed::curlLambda_(const basis& basis_,
                              const int& e) const{

    return basis_.unitTangents_.at(e); 
}

vertex Hdivmixed::curlR_(const basis& basis_,
                         const int& e1,
                         const int& e2,
                         const vertex& point) const{

    return -2 * (basis_.lambda(e1,point)*basis_.unitTangents_.at(e2) - 
                 basis_.lambda(e2,point)*basis_.unitTangents_.at(e1))/
                pow(basis_.lambda(e1,point)+basis_.lambda(e2,point),2);

}

vertex Hdivmixed::curlR_(const basis& basis_,
                         const int& e,
                         const vertex& point) const{

    return -0.5*curlR_(basis_,e,(e+2)%4,point);
}

double Hdivmixed::phie_(const basis& basis_,
                        const int& e,
                        const vertex& point) const{

    return basis_.lambda((e+1)%4,point)*basis_.lambda((e+3)%4,point)*basis_.R(e,point);

}

double Hdivmixed::phiv_(const basis& basis_,
                        const int& e,
                        const vertex& point) const{

    vertex mid1 = (basis_.corners_.at((e+3)%4) + basis_.corners_.at(e))/2.0;
    vertex mid2 = (basis_.corners_.at((e+4)%4) + basis_.corners_.at((e+1)%4))/2.0;

    return basis_.lambda((e+3)%4,point)*basis_.lambda((e+2)%4,point) - 
           basis_.lambda((e+3)%4,mid1)*basis_.lambda((e+2)%4,mid1)*
           phie_(basis_,e,point)/phie_(basis_,e,mid1) - 
           basis_.lambda((e+3)%4,mid2)*basis_.lambda((e+2)%4,mid2)*
           phie_(basis_,(e+1)%4,point)/phie_(basis_,(e+1)%4,mid2) ;

}

vertex Hdivmixed::curlphiv_(const basis& basis_,
                            const int& e,
                            const vertex& point) const{

    // Curl of nodal basis on r = 2
    // curl(DS_2)
    // phiv_i = lambda_i+2 * lambda_i+3 
    // - lambda_i+2(x_e,i)  *lambda_i+3(x_e,i)  *phi_e,i
    // - lambda_i+2(x_e,i+1)*lambda_i+3(x_e,i+1)*phi_e,i+1

    vertex mid1 = (basis_.corners_.at((e+3)%4) + basis_.corners_.at(e))/2.0;
    vertex mid2 = (basis_.corners_.at((e+4)%4) + basis_.corners_.at((e+1)%4))/2.0;

    vertex work;

    work = basis_.unitTangents_.at((e+2)%4)*basis_.lambda((e+3)%4,point) + 
           basis_.unitTangents_.at((e+3)%4)*basis_.lambda((e+2)%4,point) -
           basis_.lambda((e+2)%4,mid1)*basis_.lambda((e+3)%4,mid1)*
           phil_(basis_,e,point)/phie_(basis_,e,mid1)                          -
           basis_.lambda((e+2)%4,mid2)*basis_.lambda((e+3)%4,mid2)*
           phil_(basis_,(e+1)%4,point)/phie_(basis_,(e+1)%4,mid2) ;

    return work;
}

vertex Hdivmixed::curlphiv_star_(const basis& basis_,
                                 const int& e,
                                 const vertex& point) const{

    vertex mid1 = (basis_.corners_.at((e+3)%4) + basis_.corners_.at(e))/2.0;
    vertex mid2 = (basis_.corners_.at((e+4)%4) + basis_.corners_.at((e+1)%4))/2.0;

    return curlphiv_(basis_,e,point)/phiv_(basis_,e,basis_.corners_.at(e)) + 
           0.5*phil_(basis_,e,point)/phie_(basis_,e,mid1) + 
           0.5*phil_(basis_,(e+1)%4,point)/phie_(basis_,(e+1)%4,mid2);
}

vertex Hdivmixed::phiv_star_star_(const basis& basis_,
                                  const int& i,
                                  const vertex& point) const{

    return point - basis_.corners_.at((i+2)%4);
}

void Hdivmixed::Test(const basis& basis_){

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

            //vertex vec = phil(basis_, Case, point);
            vertex vec = phic(basis_,1,point);

            //double val = phiv_(basis_,Case, point)/phiv_(basis_,Case,basis_.corners_.at(Case));
            vertex cNormal = basis_.unitNormals_.at(Case);

            cout << cNormal[0]*vec[0] + cNormal[1]*vec[1] << " "; 

        } cout << endl<< endl;

    }
}
