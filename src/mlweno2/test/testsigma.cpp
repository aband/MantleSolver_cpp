#include "stencilpolynomial.h"
#include "tensor.h"

extern "C"{
#include "mesh.h"
#include "output.h"
}

int main(int argc, char ** argv){


    // Create pseudo mesh for testing  
    int M = 3;
    int N = 3;
    double hx = 1.0/(double)M;
    double hy = 1.0/(double)N;

    Tensor<vertex> mesh = Tensor<vertex>(2);

    mesh.setSize({M+1,N+1});

    for (int j=0; j<N+1; j++){
        for (int i=0; i<M+1; i++){
            mesh({i,j}) = {i*hx, j*hy};
        }
    }

    vector<vertex>  refcell;
    vertex cellcenter = (mesh({0,0}) + mesh({1,0}) + mesh({1,1}) + mesh({0,1}))/4;
    refcell.push_back(cellcenter);

    cellcenter = (mesh({2,0}) + mesh({3,0}) + mesh({3,1}) + mesh({2,1}))/4;
    refcell.push_back(cellcenter);

    cellcenter = (mesh({2,2}) + mesh({3,2}) + mesh({3,3}) + mesh({2,3}))/4;
    refcell.push_back(cellcenter);

    cellcenter = (mesh({0,2}) + mesh({1,2}) + mesh({1,3}) + mesh({0,3}))/4;
    refcell.push_back(cellcenter);

    double refarea = NumIntegralFace(refcell, {0,0}, {0.0,0.0}, 1.0, constFunc);

    double h = pow(refarea,0.5);

    vertex center = (mesh({0,0}) + mesh({3,0}) + mesh({0,3}) + mesh({3,3}))/4;

    // Create stencil polynomial
    vector<vector<vertex>> cornerSetSet2;
    vector<indice> order {{0,0},{1,0},{1,1},{0,1}};

    for (int j=0; j<N; j++){
    for (int i=0; i<M; i++){
        indice start = {i,j};
        indice now ; 
        vector<vertex> work;
        for (const auto& it: order){now = start + it; work.push_back(mesh({now[0],now[1]}));}
        cornerSetSet2.push_back(work);
    }}
    stencilpolynomial stencilpoly = stencilpolynomial(3,3);

    stencilpoly.setCoef(cornerSetSet2, center, h); 

    vector<double> tmp {0.037037037037043, 0.259259259259265, 0.703703703703710,
                        0.037037037037043, 0.259259259259265, 0.703703703703710,
                        0.037037037037043, 0.259259259259265, 0.703703703703710};

    Tensor<double> sol = Tensor<double>(2);

    sol.setSize({3,3});

    for (int i=0; i<9; i++) {sol(i) = tmp.at(i);}

    vertex test {0.27,0.27};

    cout << std::setprecision(10) << stencilpoly.eval(sol, test, center, h) << "  " << 0.27*0.27 << endl;

    stencilpoly.sigma(refcell, refarea, center, h);

    cout << "Current smoothness indicator is : " << stencilpoly.sigma(sol) << endl;;

    return 1;
}
