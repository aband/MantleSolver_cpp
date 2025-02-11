// manually assign everything

#include "stencilpolynomial.h"
#include "reconstruction.h"
#include "mluse.h"
#include "petsc.h"
#include "input.h"

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

    //double h = pow(refarea,0.5);
    double h = 1.0/3.0;

    vertex center = (mesh({0,0}) + mesh({3,0}) + mesh({0,3}) + mesh({3,3}))/4;

    // Create reference cell corners
    vertex add = {-h/2, -h/2};  
    refcell[0] = center + add;  
    add = {h/2, -h/2};                               
    refcell[1] = center + add;                   
    add = {h/2, h/2};                                               
    refcell[2] = center + add;                                    
    add = {-h/2, h/2};                                     
    refcell[3] = center + add;

    double refarea = NumIntegralFace(refcell, {0,0}, {0.0,0.0}, 1.0, constFunc);

    return 1;
}
