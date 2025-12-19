#ifndef CLEAN_SOLVER_H_
#define CLEAN_SOLVER_H_

#include <petsc.h>
#include <iostream>
#include <ctime>
#include <chrono>
#include <map>
#include "integral.h"
#include "eutectic_rescaled.h"
//#include "input.h"
#include "util.h"

// MFEM parameter header file
#include "myFunc.h"

#include "basis.h"
#include "Hdivmixed.h"
#include "brmixed.h"
#include "shape.h"

// =========================================================================
// The boundary value data structure contains
// 1. global index of degree of freedom and global index of element
// 2. a pair object pairing local degree of freedom and value
struct bndryInfo{
    int        localDOF;
    indice     globalElem;

    double     essenval;
    double     naturval;
};

// Map global indice with boundary values
using bndryVal = std::map<int, bndryInfo>;

using bndryValGroup = std::map<std::string, std::vector<bndryVal>>;

typedef struct{

    std::vector<double> A;
    std::vector<double> B;
    double C;
    std::vector<double> f;

} LocMat;

typedef struct{

    std::vector<double> edgeporo;
    std::vector<double> cellporo;
    double aveporo;

} poroSet;

typedef struct{

    Mat M, Kg, B, Bg, C;
    Vec g, source, neum;

    // Later added vectors
    Vec F, G;
    Vec x, y;

} ReducedSys;

bndryType bMarker(const MeshInfo& mi, std::vector<int> work,
                  const std::string& name,
                  const std::vector<double>& parameter);

template <typename T>
int CreateRefMap(T& funcSp, 
                 const MeshInfo& mi, 
                 int * refArray,
                 std::map<int, int>& refMapNatur,
                 int& EssenDOFCount,
                 int& NaturDOFCount,
                 const std::vector<double>& parameter){

    int essenCount = 0; 
    int naturCount = 0;
    int interCount = 0;

    // ! loop through all dofs 
    // Natural boundary dofs are different from Essential boundary dofs
    // For the fact that natrual boundary dofs participate in the left hand side matrix
    // and also right hand side vector
    // Hence it requires two different index system.

    for (int dof=0; dof<funcSp.getDOF(); dof++){

        std::vector<int> work = funcSp.GlobalToLocalMapBndry(mi, dof);

        bndryType bt = bMarker(mi, work, funcSp.name, parameter);

        if (bt == dirichlet){
            refArray[dof] = essenCount;
            essenCount ++;
        } else if (bt == neumann){
            refArray[dof] = interCount;
            interCount ++;
            refMapNatur.insert(std::make_pair(dof, naturCount));
            naturCount ++;
        } else {
            // bt == missed
            // Currently only mixing essential and natural bondary conditions
            // Hence, essential and natural boundary should consists of all the boundarys
            refArray[dof] = interCount;
            interCount ++;
        }

    }

    EssenDOFCount = essenCount;
    NaturDOFCount = naturCount;

    return 1;
}

template <typename T>
int AssignLocRedSys(ReducedSys& redsys,
                    LocMat& loc,
                    int * ref,
                    const map<int, int>& refNatur,
                    const MeshInfo& mi,
                    const bndryVal& bndryAll,
                    const indice& global,
                    T& funcSp,
                    const std::vector<double>& parameter){

    // ! Get global index of local dofs 
    //const std::vector<int> elemDofs = funcSp.LocalToGlobal(mi, global);
    const std::vector<int> elemDofs = funcSp.LocalGlobalMap(mi, global);

    const int idxn = FlatIndic(mi, global);

    for (int row=0; row<elemDofs.size(); row++){
        // View it as the row index
        const int idxm = ref[elemDofs.at(row)];
        const double valB = loc.B.at(row);

//        if (funcSp.onBndry(mi, elemDofs.at(row))){
        if (bMarker(mi,funcSp.GlobalToLocalMapBndry(mi,elemDofs.at(row)),funcSp.Name(),parameter) == dirichlet){

            // This dof is a dirichlet dof on the boundary
            // Should be assign to Bg
            PetscCall(MatSetValues(redsys.Bg, 1, &idxm, 1, &idxn, &valB, 
                                   ADD_VALUES));

            // At the same time insert essential boundary value to rhs vector
            auto itFind = bndryAll.find(elemDofs.at(row));
            if (itFind != bndryAll.end()){
                const bndryInfo& tmp = bndryAll.at(elemDofs.at(row));
                PetscCall(VecSetValues(redsys.g, 1, &idxm, &tmp.essenval, INSERT_VALUES));
            }
        } else {
            // This dof is not on the boundary
            // Should be assigned to B instead
            PetscCall(MatSetValues(redsys.B , 1, &idxm, 1, &idxn, &valB, 
                                   ADD_VALUES));

            // This dof is not on the boundary
            // This dof will contribute to source term
            double vals = loc.f.at(row);
            PetscCall(VecSetValues(redsys.source, 1, &idxm, &vals, ADD_VALUES));

            for (int col=0; col<elemDofs.size(); col++){
                const int cidxn  = ref[elemDofs.at(col)];
                const double val = loc.A.at(row+col*elemDofs.size()); 

//                if (funcSp.onBndry(mi,elemDofs.at(col))){
                  if (bMarker(mi,funcSp.GlobalToLocalMapBndry(mi,elemDofs.at(col)),funcSp.Name(),parameter) == dirichlet){

                    // It is a non bndry dof - bndry dof interaction
                    // Val assigned to M
                    PetscCall(MatSetValues(redsys.Kg, 1, &idxm, 1, &cidxn, 
                                           &val, ADD_VALUES));
                } else {
                    // It is a non bndry dof - non bndry dof interaction
                    // Val assigned to M
                    PetscCall(MatSetValues(redsys.M, 1, &idxm, 1, &cidxn, 
                                           &val, ADD_VALUES));
                }
            }

            // Add correction to natural dof ============================================================================
            auto itFind = refNatur.find(elemDofs.at(row));

            if (itFind != refNatur.end()){
                const bndryInfo& tmp = bndryAll.at(elemDofs.at(row));
                PetscCall(VecSetValues(redsys.neum, 1, &idxm, &tmp.naturval, INSERT_VALUES));
            }

            // ==========================================================================================================
        }
    }

    return 1;
}

template <typename T>
vector<vertex> ExtractVelocity(Vec * sol, Vec * g,
                               int *refmap,
                               const MeshInfo& mi,
                               vector<vertex> points,
                               const indice& gCell,
                               T& funcSp,
                               basis& mybasis,
                               const std::vector<double>& parameter){

    std::vector<vertex> work;
    work.resize(points.size());

    //PetscScalar *valuesSol;
    //PetscScalar *valuesg;

    double *valuesSol;
    double *valuesg;

    VecGetArray(*sol, &valuesSol);
    VecGetArray(*g, &valuesg);

    mybasis.GetCorners(mi, gCell);

    // !Get global indiex of the local dofs in specific and correct order
    const std::vector<int> elemDofs = funcSp.LocalGlobalMap(mi, gCell);

    for (int g=0; g<points.size(); g++){

        // Initialize interpolated value
        work.at(g) = {0.0,0.0};

        std::vector<vertex> basisVal = funcSp.EvaluateAll(mybasis, points.at(g));

        // Reconstruction of value with element basis
        for (int k=0; k<elemDofs.size(); k++){

            if (bMarker(mi,funcSp.GlobalToLocalMapBndry(mi,elemDofs.at(k)),funcSp.name, parameter) == dirichlet){
                work.at(g) += valuesg[refmap[elemDofs.at(k)]] * basisVal.at(k);
            } else {
                work.at(g) += valuesSol[refmap[elemDofs.at(k)]] * basisVal.at(k);
            }
        }
    }

    VecRestoreArray(*sol, &valuesSol);
    VecRestoreArray(*g, &valuesg);

    return work;
}

// Extract velocity on given edge gauss points set
// Vertical and Horizontal edges
template <typename T>
int ExtractVelocityEdge(vector<vertex>& edgeVelocity,
                        const vector<vertex>& edgegaussp,
                        const MeshInfo& mi,
                        int* refmap, 
                        Vec * sol, Vec * g, 
                        T& funcSp, basis& mybasis){

    indice gCell;

    edgeEnds<vertex> edgeEndsVertex;
    edgeEnds<indice> edgeEndsIndice;

    vector<vertex> gaussp;
    gaussp.resize(3);

    vector<vertex> velgauss;
    velgauss.resize(3);

    int M = mi.MPIglobalCellSize[0];
    int N = mi.MPIglobalCellSize[1];

    int tolvert = (M+1)*N*3;
    int tolhori = M*(N+1)*3;

    int toledgegauss = tolvert + tolhori;

    // Vertical points first
    for (int j=0; j<N  ; j++){
    for (int i=0; i<M+1; i++){

        int dof = (j*(M+1) + i)*3;

        velgauss.clear();
        velgauss.resize(3);
        gaussp.clear();
        gaussp.resize(3);

        if (i == M){
            // Right boundary
            gCell = {i-1,j};
        } else {
            gCell = {i,j};
        }

        for (int g=0; g<3; g++){
            gaussp.at(g) = edgegaussp.at(dof+g);
        }

        velgauss = ExtractVelocity(sol, g, refmap, mi, gaussp, gCell, funcSp, mybasis, {1});

        for (int g=0; g<3; g++){
            edgeVelocity.at(dof+g) = velgauss.at(g);
        }

    }}

    for (int j=0; j<N+1; j++){
    for (int i=0; i<M  ; i++){

        int dof = tolvert + (j*M+i)*3;

        velgauss.clear();
        velgauss.resize(3);
        gaussp.clear();
        gaussp.resize(3);

        if (j == N){
            // Right boundary
            gCell = {i,j-1};
        } else {
            gCell = {i,j};
        }

        for (int g=0; g<3; g++){
            gaussp.at(g) = edgegaussp.at(dof+g);
        }

        velgauss = ExtractVelocity(sol, g, refmap, mi, gaussp, gCell, funcSp, mybasis, {1});

        for (int g=0; g<3; g++){
            edgeVelocity.at(dof+g) = velgauss.at(g);
        }

    }}

    return 1;
}

class DarcyStokes{
    public:
        DarcyStokes(const MeshInfo& mi, PhysProperty * pp, const std::vector<double>& param) {init(mi, pp, param);};
        ~DarcyStokes() 
         {delete refArrayStokesEssen_; 
          delete refArrayDarcyEssen_;};

        vector<vertex> StokesVel; 
        vector<vertex> DarcyVel;

        int init(const MeshInfo& mi,
                 PhysProperty * pp,
                 const std::vector<double>& param);

        int Assemble(const MeshInfo& mi,
                     const vector<double>& edgeporo,
                     const vector<double>& cellporo,
                     const vector<double>& averporo,
                     double theta,
                     PhysProperty * pp);

        int CreateCoupledSystem();

        int Solve(int maxIter, double tolUzawa);

        int ReconstructEdgeVel(const vector<vertex>& edgegaussp,
                               const MeshInfo& mi);

        // Printing functions
        int printBndryAll();

        int showMatrix();

        // Some supplemental set for local dofs
        // normal dof 
        std::set<int> top_normal_stokes {11,7,6};
        std::set<int> left_normal_stokes {0,3,8};
        std::set<int> right_normal_stokes {1,2,10};
        std::set<int> bottom_normal_stokes {4,5,9};

        std::vector<std::set<int>> normal_stokes {left_normal_stokes,
                                                  bottom_normal_stokes,
                                                  right_normal_stokes,
                                                  top_normal_stokes};
        // tangent dof
        std::set<int> top_tang_stokes {3,2};
        std::set<int> left_tang_stokes {7,4};
        std::set<int> right_tang_stokes {5,6};
        std::set<int> bottom_tang_stokes {0,1};

        std::vector<std::set<int>> tang_stokes {left_tang_stokes,
                                                bottom_tang_stokes,
                                                right_tang_stokes,
                                                top_tang_stokes};

    private:

        // Total number of elements
        int totalElem;

        int ComputeEssenBndryAll(const MeshInfo& mi,
                                 PhysProperty * pp,
                                 const std::vector<double>& param);

        int ComputeNaturBndryAll(const MeshInfo& mi,
                                 PhysProperty * pp,
                                 const std::vector<double>& param);

        double AssignBndrySupVal(const vertexSet& edgeCorner,
                                 const vertex& nu,
                                 PhysProperty * pp);

        std::array<double, 2> AssignBndryValsDarcy(const vertexSet& edgeCorner,
                                                   const vertex& nu,
                                                   const double& len,
                                                   const int& edge,
                                                   PhysProperty * pp);

        int computeEssenVals(const MeshInfo& mi, int i, int j, int edge, PhysProperty * pp);
        int computeNaturVals(const MeshInfo& mi, int i, int j, int edge, PhysProperty * pp);

        int AssignLocMatDarcy(const MeshInfo& mi,  LocMat& loc, double theta, const poroSet& poro, PhysProperty * pp);
        int AssignLocMatStokes(const MeshInfo& mi, LocMat& loc, double theta, const poroSet& poro, PhysProperty * pp);
        int AssignLocMatCouple(const MeshInfo& mi, double& k,   double theta, const poroSet& poro, PhysProperty * pp);

        // Assembler for the interior cells
        int AssignLocRedSysDarcy(LocMat& loc,
                                 int * ref,
                                 const MeshInfo& mi,
                                 const indice& global);

        int AssignLocRedSysStokes(LocMat& loc,
                                  int * ref,
                                  const MeshInfo& mi,
                                  const indice& global);

        // Assembler for the boundary cells
//        int AssignLocRedSysDarcy(LocMat& loc,
//                                 int * ref,
//                                 const MeshInfo& mi,
//                                 const bndryVal& bndryAll,
//                                 const indice& global);

//        int AssignLocRedSysStokes(LocMat& loc,
//                                  int * ref,
//                                  const MeshInfo& mi,
//                                  const bndryVal& bndryAll,
//                                  const indice& global);

        int PrepareReducedSys(ReducedSys& redsys, 
                              int reducedDOF, int bndrySize, 
                              int Adnz, int Aonz, int Bdnz, int Bonz);

        // Function basis
        basis     basis_;
        Hdivmixed hdiv_;
        BRMixed   br_;

        /* ====================================================  
         Boundary values containing four edges numberred 0-3
         The numbering order is the same as the element edge number
		   0 -- left
		   1 -- bottom
		   2 -- right
		   3 -- top
			The corner dofs are defined and dealt with separately
		  ==================================================== */ 

        bndryVal bndryStokesAll;
        bndryVal bndryDarcyAll;

        /**!
         * Reduced linear system excluding essential boundary conditions
         */
        ReducedSys reducedDarcy_;
        ReducedSys reducedStokes_;

        // Struct holding computed result
        ReducedSys result;

        /**!
         * Coupling matrix.
         */
        Mat K;

        /**!
         * Create boundary dof reference mapping
         */
        int bndryDOFStokesEssen_ = 0.0;
        int bndryDOFDarcyEssen_  = 0.0;
 
        int bndryDOFStokesNatur_ = 0.0;
        int bndryDOFDarcyNatur_ = 0.0;

        /**!
         * Reference map tell what kind of boundary condition dof belongs to
         */
        int * refArrayStokesEssen_;
        int * refArrayDarcyEssen_;
 
        map<int,int> refArrayStokesNatur_;
        map<int,int> refArrayDarcyNatur_;

        // Size of the remaining dofs
        int reducedDOFStokes;
        int reducedDOFDarcy;
};

#endif
