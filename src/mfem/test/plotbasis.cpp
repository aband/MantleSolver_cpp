// Test basis soely
#include <iostream>
#include <petsc.h>
#include "integral.h"
#include "input.h"
#include "Hdivmixed.h"
#include "brmixed.h"
#include "assemble.h"
#include "util.h"
#include "bndry.h"
#include "solve.h"
#include "error.h"

extern "C"{
#include "mesh.h"
#include "output.h"
}

using namespace std;

vertex bilinearMap(const vertexSet& corners, const vertex& ref){

    vertex target {0,0};

    target[0] = (1-ref[0])*(1-ref[1])*corners.at(0)[0] + ref[0]*(1-ref[1])*corners.at(1)[0] + ref[0]*ref[1]*corners.at(2)[0] + (1-ref[0])*ref[1]*corners.at(3)[0];
    target[1] = (1-ref[0])*(1-ref[1])*corners.at(0)[1] + ref[0]*(1-ref[1])*corners.at(1)[1] + ref[0]*ref[1]*corners.at(2)[1] + (1-ref[0])*ref[1]*corners.at(3)[1];

    return target;
}

int main(int argc, char ** argv){

    // Declare classes for shape functions
    basis * basis_   = new basis();
    Hdivmixed * hdiv = new Hdivmixed();
    BRMixed * br     = new BRMixed();

    int dof=0;
    PetscOptionsGetInt(NULL,NULL,"-dof",&dof,NULL);

    int seed = 21;
    double h = 1.0/(double) (seed-1) ;

    // Define a mapping from reference [1 0, 0 1] to a quad
    vertex v0 {0.1,-0.2};
    vertex v1 {0.8,0.1};
    vertex v2 {1.2,0.95};
    vertex v3 {-0.05, 1.03};

    vertexSet corners = {v0,v1,v2,v3};

    // reference element
    //basis_->GetCorners(mi, {0,0});

    basis_->GetCorners(corners);
 
    FILE *fx = fopen("gridX.dat","w");
    FILE *fy = fopen("gridY.dat","w");
    FILE *val1 = fopen("R1.dat","w");
    FILE *val2 = fopen("R.dat","w");
    FILE *val3 = fopen("phie.dat","w");
    FILE *val4 = fopen("phiv.dat","w");

    for (int j=0; j<seed; j++){
    for (int i=0; i<seed; i++){
        vertex sample {i*h, j*h};
        vertex target = bilinearMap(corners, sample);;

        fprintf(fx, "%f ", target[0]);
        fprintf(fy, "%f ", target[1]);
        fprintf(val1, "%f ", basis_->R(0,2,target));
        fprintf(val2, "%f ", basis_->R(0,target));
        fprintf(val3, "%f ", hdiv->phie(*basis_,0,target)/hdiv->phie(*basis_,0,(v3+v0)/2));
		  fprintf(val4, "%f ", hdiv->phiv(*basis_,0,target)/hdiv->phiv(*basis_,0,(v0)));
    }}

    fclose(fx);
    fclose(fy);
    fclose(val1);
    fclose(val2);
    fclose(val3);
	 return 0;
}
