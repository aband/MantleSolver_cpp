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
cout << h << endl;
    stencilpoly.setCoef(cornerSetSet2, center, h); 

    stencilpoly.printCoef();

    vector<double> tmp {0.037037037037043, 0.259259259259265, 0.703703703703710,
                        0.037037037037043, 0.259259259259265, 0.703703703703710,
                        0.037037037037043, 0.259259259259265, 0.703703703703710};

    Tensor<double> sol = Tensor<double>(2);

    sol.setSize({3,3});

    for (int i=0; i<9; i++) {sol(i) = tmp.at(i);}

    vertex test {0.27,0.27};

    polynomial collapse = stencilpoly.createCollapsePoly(sol);

    cout << std::setprecision(10) << stencilpoly.eval(sol, test, center, h) << "  " << 0.27*0.27 << endl;

    cout << "Collapse polynomial result : " << collapse.eval((test[0]-center[0])/h, (test[1]-center[1])/h) << endl;

    stencilpoly.sigma(refcell, refarea, center, h);

    cout << "Current smoothness indicator is : " << stencilpoly.sigma(sol) << endl;;

    cout << "Smoothness indicator calculated using collapse is : " << stencilpoly.sigma(collapse, refcell, refarea, center, h) << endl;

    // Test part 2
    Tensor<double> sol1 = Tensor<double>(2);
    sol1.setSize({2,2});
    Tensor<double> sol2 = Tensor<double>(2);
    sol2.setSize({2,2});
    Tensor<double> sol3 = Tensor<double>(2);
    sol3.setSize({2,2});
    Tensor<double> sol4 = Tensor<double>(2);
    sol4.setSize({2,2});
 
    sol1({0,0}) = sol({0,0});
    sol1({1,0}) = sol({1,0});
    sol1({1,1}) = sol({1,1});
    sol1({0,1}) = sol({0,1});

    sol2({0,0}) = sol({1,0});
    sol2({1,0}) = sol({2,0});
    sol2({1,1}) = sol({2,1});
    sol2({0,1}) = sol({1,1});

    sol3({0,0}) = sol({1,1});
    sol3({1,0}) = sol({2,1});
    sol3({1,1}) = sol({2,2});
    sol3({0,1}) = sol({1,2});

    sol4({0,0}) = sol({0,1});
    sol4({1,0}) = sol({1,1});
    sol4({1,1}) = sol({1,2});
    sol4({0,1}) = sol({0,2});

    stencilpolynomial stencilpoly1 = stencilpolynomial(2,2);
    stencilpolynomial stencilpoly2 = stencilpolynomial(2,2);
    stencilpolynomial stencilpoly3 = stencilpolynomial(2,2);
    stencilpolynomial stencilpoly4 = stencilpolynomial(2,2);

    vertex center_test2 {0.5,0.5};
    vector<vertex> refcell_test2;
    refcell_test2.push_back(mesh({1,1}));
    refcell_test2.push_back(mesh({2,1}));
    refcell_test2.push_back(mesh({2,2}));
    refcell_test2.push_back(mesh({1,2}));

    vector<vector<vertex>> cornerSet_test2_1;
    vector<vector<vertex>> cornerSet_test2_2;
    vector<vector<vertex>> cornerSet_test2_3;
    vector<vector<vertex>> cornerSet_test2_4;

    for (int j=0; j<2; j++){
    for (int i=0; i<2; i++){
        indice start = {i,j};
        indice now ; 
        vector<vertex> work1;
        vector<vertex> work2;
        vector<vertex> work3;
        vector<vertex> work4;
        for (const auto& it: order){
             now = start + it; 
             work1.push_back(mesh({now[0],now[1]}));
             work2.push_back(mesh({now[0]+1,now[1]}));
             work3.push_back(mesh({now[0]+1,now[1]+1}));
             work4.push_back(mesh({now[0],now[1]+1}));
        }
        cornerSet_test2_1.push_back(work1);
        cornerSet_test2_2.push_back(work2);
        cornerSet_test2_3.push_back(work3);
        cornerSet_test2_4.push_back(work4);

    }}
 
    stencilpoly1.setCoef(cornerSet_test2_1, center_test2, h); 
    stencilpoly2.setCoef(cornerSet_test2_2, center_test2, h); 
    stencilpoly3.setCoef(cornerSet_test2_3, center_test2, h); 
    stencilpoly4.setCoef(cornerSet_test2_3, center_test2, h); 

    stencilpoly1.sigma(refcell_test2, h*h, center_test2, h);
    stencilpoly2.sigma(refcell_test2, h*h, center_test2, h);
    stencilpoly3.sigma(refcell_test2, h*h, center_test2, h);
    stencilpoly4.sigma(refcell_test2, h*h, center_test2, h);

    cout << stencilpoly1.sigma(sol1) << endl;
    cout << stencilpoly2.sigma(sol2) << endl;
    cout << stencilpoly3.sigma(sol3) << endl;
    cout << stencilpoly4.sigma(sol4) << endl;

    return 1;
}
