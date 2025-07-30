#include <iostream> 
#include <cmath>
#include <cstdlib>

#include "../polynomial.h"

#include <chrono>

#include "reconstruction.h"
#include "mluse.h"
#include "petsc.h"
#include "input.h"

extern "C"{
#include "mesh.h"
#include "output.h"
}

// Horner's method
static double horner(double x, const double* coef, int degree) {
  if(abs(x) <= 1) {

    double val = coef[degree];
    for(int i = degree-1; i >= 0; i--) {
      val = val*x + coef[i];
    }
    return val;

  } else {
   
    double val = coef[0];
    for(int i = 1; i <= degree; i++) {
      val = val/x + coef[i];
    }
    return pow(x,degree) * val;
    
  }
}

// 1D Horner's method for polynomial derivative evaluation (all up to der) for p(x/h)
static void horner_der(int der, double* val, double x, double h, 
					        const double* coef, int degree) {
  double xx = x/h;

//  for (int t=0; t<degree+1; t++){
//  cout << coef[t] << "  " ;
//  }
//cout << endl;
  for(int i = 0; i<= der; i++) val[i] = 0;
  
  for(int i = degree; i >= 0; i--) {
    for(int d=der; d>=1; d--) val[d] = val[d]*xx + d*val[d-1];
    val[0] = val[0]*xx + coef[i];
  }

  for(int d=1; d<=der; d++) val[d] /= pow(h,d);
}

// 1D Horner's method for all polynomial derivative evaluation for p(x/h)
static void horner_der(double* val, double x, double h, 
					        const double* coef, int degree) {
  double xx = x/h;
  
  for(int i = 0; i<= degree; i++) val[i] = 0;
  
  for(int i = degree; i >= 0; i--) {
    for(int d=degree - i; d>=1; d--) val[d] = val[d]*xx + d*val[d-1];
    val[0] = val[0]*xx + coef[i];
  }

  for(int d=1; d<=degree; d++) val[d] /= pow(h,d);
}

// Evaluate val[m + (derX+1)*n] = D_x^m D_y^n p(x), p(x) = Sum_ij c_ij (x-x0)^i (y-y0)^j / h^(i+j)
void polynomial2D_ders(int derX, int derY, double* val,
			              double x, double y, double x0, double y0, double h,
			              int polyn_degree, double* my_coef) {
  double xx0 = x-x0;
  double yy0 = y-y0;

  double valX[derX+1];
  double valY[derY+1];
  double yCoef[derX+1][polyn_degree+1];
  
  // Horner's method in x, for each power of y
  int sz = polyn_degree+1;
  int start = 0;
  for(int j = 0; j <= polyn_degree; j++) {
    horner_der(derX,valX,xx0,h,&my_coef[start],sz-1);
    for(int d = 0; d <= derX; d++) 
	 {yCoef[d][j] = valX[d];}
    start += sz;
    //sz--;
  }

  // Horner's method in y
  for(int dX = 0; dX <= derX; dX++) {
    horner_der(derY,valY,yy0,h,yCoef[dX],polyn_degree);
    for(int dY = 0; dY <= derY; dY++) {
      val[dX + (derX+1)*dY] = valY[dY];
    }
  }
}

int main(int argc, char ** argv){

    // A n*n tensor product polynomial
    double my_coef[5*5];

    double val[5*5];

    for (int i=0; i<25; i++){
        my_coef[i] = (double)i;
        val[i]     = 0.0;
    }

/*
    my_coef[9] = 0.0;
    my_coef[13] = 0.0;
	 my_coef[14] = 0.0;
	 my_coef[17] = 0.0;
	 my_coef[18] = 0.0;
	 my_coef[19] = 0.0;
	 my_coef[21] = 0.0;
	 my_coef[22] = 0.0;
	 my_coef[23] = 0.0;
	 my_coef[24] = 0.0;
*/

	 auto start = std::chrono::steady_clock::now();
	 for (int it = 0; it < 1000000; it++){
    polynomial2D_ders(4,4,val,1.0,1.0,0.0,0.0,0.2,4,my_coef);}
	 auto end = std::chrono::steady_clock::now();
	 auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end-start);
    cout << "Benchmark Test: " << duration.count() << " ms." << endl;

    // Benchmark code
    for (int i=0; i<25; i++){
        std::cout << val[i] << " " ;
		  if (i%5 == 4){
            std::cout << std::endl;
		  }
    }

    cout << endl; 

    // =================================
    polynomial testp = polynomial(5,5);
    testp.setCoef(my_coef);
    //testp.printCoef();

    Tensor<double> der;
    der.setSize({5,5});

	 start = std::chrono::steady_clock::now();
	 for (int it = 0; it < 1000000; it++){
    testp.evalDer(4,4,1.0,1.0,0.2,der);
    }
	 end = std::chrono::steady_clock::now();
	 duration = std::chrono::duration_cast<std::chrono::milliseconds>(end-start);
    cout << "My function Test: " << duration.count() << " ms." << endl;

    for (int i=0; i<25; i++){
        std::cout << der(i) << " " ;
		  if (i%5 == 4){
            std::cout << std::endl;
		  }
    }

    // ==================================
	  
    cout << endl << "Test of stencil sigmas ... " << endl;

    // Initializing petsc function
    PetscErrorCode ierr;
    PetscMPIInt   size,rank;
    PetscInitialize(&argc, &argv, NULL, NULL);

    MPI_Init(NULL,NULL);
    MPI_Comm_size(PETSC_COMM_WORLD,&size);
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    int M = 240, N = 80;
    ierr = PetscOptionsGetInt(NULL,NULL,"-M",&M,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL,NULL,"-N",&N,NULL);CHKERRQ(ierr);
    double L = 3, H = 1;
    double xstart = 0.0, ystart = 0.0;
    PetscCall(PetscOptionsGetReal(NULL,NULL,"-L",&L,NULL));
    PetscCall(PetscOptionsGetReal(NULL,NULL,"-H",&H,NULL));
    PetscCall(PetscOptionsGetReal(NULL,NULL,"-xstart", &xstart, NULL));
    PetscCall(PetscOptionsGetReal(NULL,NULL,"-ystart", &ystart, NULL));

    int stencilWidthMesh = 5; // Ghost layer thickness for vertex
    int stencilWidthU = 3;    // Ghost layer thickness for cell

    DM dmu;
    DM dmMesh;

    int meshType = 0; 
    PetscCall(PetscOptionsGetInt(NULL,NULL,"-meshtype",&meshType,NULL));

    double dscale = 1.0;
    PetscCall(PetscOptionsGetReal(NULL,NULL,"-scale",&dscale,NULL));

    L/=dscale;
    H/=dscale;

    // Create dmMesh
    PetscCall(DMDACreate2d(PETSC_COMM_WORLD, 
    DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DMDA_STENCIL_BOX, 
    M, N, PETSC_DECIDE, PETSC_DECIDE, 2, stencilWidthMesh, NULL, NULL, &dmMesh));
    PetscCall(DMSetFromOptions(dmMesh));              
    PetscCall(DMSetUp(dmMesh));

    // Create dmU
    PetscCall(DMDACreate2d(PETSC_COMM_WORLD, 
    DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DMDA_STENCIL_BOX, 
    M, N, PETSC_DECIDE, PETSC_DECIDE, 1, 
    stencilWidthU, NULL, NULL, &dmu));
    PetscCall(DMSetFromOptions(dmu));              
    PetscCall(DMSetUp(dmu));     

	 MeshParam mp; 
    mp.xstart = xstart;
    mp.ystart = ystart;
    mp.L = L;
    mp.H = H;

    Vec globalmesh;
    // Create global vector containing mesh
    PetscCall(DMCreateGlobalVector(dmMesh, &globalmesh));

    switch(meshType){
        case 0: CreateFullMesh(dmMesh, &globalmesh, &mp); break;
        case 1: LogicRectMesh(dmMesh, &globalmesh, &mp);  break;
        case 2: RefineMesh(dmMesh, &globalmesh, &mp);
        //case 2: TestControlMeshSecond(dmCell,L,H); break;
        //case 3: TestControlMeshThird(dmCell,L,H);  break;
    }

    MeshInfo mi;

    ReadMeshPortion(dmMesh, &globalmesh, mi.lmesh);

    AssignValuesMeshInfo(mi, dmMesh, dmu);

    mi.L = L;
    mi.H = H;

    double h0 = sqrt((L*H)/(double)(M*N));

    // =====================================================================
    multilevel ml = multilevel();

    cout << "Efficiency Test on grid : " << M << "  " << N << endl;

	 start = std::chrono::steady_clock::now();
	 ml.addLevel("(5,5)",{5,5},mi);
	 end = std::chrono::steady_clock::now();
	 duration = std::chrono::duration_cast<std::chrono::milliseconds>(end-start);
    cout << "(5,5) level created. Using " << duration.count() << " ms." << endl;

    // Pick a stencil
    int startM = M/2-2;
    int startN = N/2-2;

    // Create a stencil polynomial
    stencilpolynomial stenp = stencilpolynomial(5);

    vector<vector<vertex>> cornerSet;
    vector<vertex> refcell;
    vertex center;

    // Extract a 5*5 stencil
    for (int j=0; j<5; j++){
        for (int i=0; i<5; i++){
            indice global {i+startM, 
                           j+startN};
            vector<vertex> cellCornerSet = extractCorners(mi, global);
            cornerSet.push_back(cellCornerSet);
        } 
    }

    return 1;
}
