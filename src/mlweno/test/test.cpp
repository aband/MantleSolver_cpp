#include "stencilpolynomial.h"
#include "tensor.h"

extern "C"{
#include "mesh.h"
#include "output.h"
}

int main(int argc, char ** argv){

    polynomial * t2 = new polynomial(2,1);
    double coef2[2] = {1,2};

    t2->setCoef(coef2);

    cout << t2->eval(1,0) << endl;
    cout << t2->eval(0,1) << endl;

    double der = 1;

    double work[2] = {0.0,0.0};

    computeDerivative(der,2,1,1,coef2,work);

    //for (int i=0; i<2; i++){
    //    cout << work[i] << " " ;
    //} cout << endl;

    // =======================================================

    int degreex = 2; 
    int degreey = 3;

    polynomial * testpoly = new polynomial(degreex,degreey);

    double testcoef[6] = {1,2,3,4,5,6};

    testpoly->setCoef(testcoef);
    //testpoly->printCoef();

    //cout << testpoly->eval(1,0) << endl;
    //cout << testpoly->eval(0,1) << endl;

    int derx = 1;
    int dery = 2;

    Tensor<double> derTensor = Tensor<double>(2);
    derTensor.setSize({derx+1, dery+1});

    cout << endl;

    testpoly->evalDer(derx,dery, .1,.1,1.0,derTensor);

    cout << endl;

    for (int dy=0; dy<dery+1; dy++){
        for (int dx=0; dx<derx+1; dx++){
            cout << derTensor({dx,dy}) << " " ;

        }cout << endl;
    }

    // =======================================================
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

    vector<vertex> cornerSet;
    cornerSet.push_back(mesh({1,1}));
    //cout << mesh({1,1})[0] << "  " << mesh({1,1})[1] << endl;
    cornerSet.push_back(mesh({2,1}));
    //cout << mesh({2,1})[0] << "  " << mesh({2,1})[1] << endl;
    cornerSet.push_back(mesh({2,2}));
    //cout << mesh({2,2})[0] << "  " << mesh({2,2})[1] << endl;
    cornerSet.push_back(mesh({1,2}));
    //cout << mesh({1,2})[0] << "  " << mesh({1,2})[1] << endl;

    vector<vertex> cornerSet1;
    cornerSet1.push_back(mesh({1,0}));
    cornerSet1.push_back(mesh({2,0}));
    cornerSet1.push_back(mesh({2,1}));
    cornerSet1.push_back(mesh({1,1}));

    vector<vertex> cornerSet2;
    cornerSet2.push_back(mesh({2,0}));
    cornerSet2.push_back(mesh({3,0}));
    cornerSet2.push_back(mesh({3,1}));
    cornerSet2.push_back(mesh({2,1}));

    vector<vertex> cornerSet3;
    cornerSet3.push_back(mesh({2,1}));
    cornerSet3.push_back(mesh({3,1}));
    cornerSet3.push_back(mesh({3,2}));
    cornerSet3.push_back(mesh({2,2}));

    vertex center = (mesh({1,0}) + mesh({3,0}) + mesh({3,2}) + mesh({1,2})) / 4.0;

    degreex = 2;
    degreey = 2;

    polynomial testint = polynomial(degreex,degreey);

    double intcoef[4] = {0,0,0,1};

    testint.setCoef(intcoef);

    cout << cornerSet.size() << endl;

    cout <<std::setprecision(20) << "poly integral : "  << polyNumIntegralFace(cornerSet, hx, center, testint) << endl; 
    cout <<std::setprecision(20) << "Previously defined integral : " << NumIntegralFace(cornerSet, {1,1}, center, hx, basePoly) << endl;

    cout << endl;
    // poly num matched with previously defined function
    stencilpolynomial teststencilpoly = stencilpolynomial(2,2);

    vector<vector<vertex>> cornerSetSet;
    cornerSetSet.push_back(cornerSet);
    cornerSetSet.push_back(cornerSet1);
    cornerSetSet.push_back(cornerSet2);
    cornerSetSet.push_back(cornerSet3);

    teststencilpoly.setCoef(cornerSetSet ,center, hx);

    teststencilpoly.printCoef();

    // ===========================================================================================================================
    cout << endl; 
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
    stencilpolynomial teststencilpoly2 = stencilpolynomial(3,3);

    center = (mesh({0,0}) + mesh({3,0}) + mesh({0,3}) + mesh({3,3}))/4;
    double h = hx;

    center = {0.1,0.123};

    h = 0.1;

    teststencilpoly2.setCoef(cornerSetSet2 ,center, h);
    teststencilpoly2.printCoef();

    vector<double> tmp {0.037037037037043, 0.259259259259265, 0.703703703703710,0.037037037037043, 0.259259259259265, 0.703703703703710,0.037037037037043, 0.259259259259265, 0.703703703703710};

    Tensor<double> sol = Tensor<double>(2);

    sol.setSize({3,3});

    for (int i=0; i<9; i++) {sol(i) = tmp.at(i);}

    vertex test {0.27,0.27};

    cout << std::setprecision(15) << teststencilpoly2.eval(sol, test, center, h) << "  " << 0.27*0.27 << endl;

    // Test for smoothness indicator
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

    // ===============================================================================
    Tensor<double> testt1 = Tensor<double>(2);
    Tensor<double> testt2 = Tensor<double>(2);
    Tensor<double> testt3 = Tensor<double>(2);

    testt1.setSize({2,2});
    testt2.setSize({2,2});
    testt3.setSize({2,2});

    for (int i=0; i<testt1.getSize(); i++){
        testt1(i) = i;
        testt2(i) = 3;
    }

    Tensor_add(testt1,testt2,testt3);
    cout << endl;
    for (int i=0; i<testt1.getSize(); i++){
        cout << testt3(i) << endl;
    }

    Tensor_multi(testt1,testt2,testt3);
    cout << endl;
    for (int i=0; i<testt1.getSize(); i++){
        cout << testt3(i) << endl;
    }



    return 0;
}
