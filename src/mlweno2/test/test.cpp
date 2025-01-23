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

    int derx = 2;
    int dery = 3;

    Tensor<double> derTensor = Tensor<double>(2);
    derTensor.setSize({derx+1, dery+1});

    cout << endl;

    testpoly->evalDer(derx,dery,degreex, degreey, .1,.1,1.0,derTensor);

    cout << endl;

    for (int dy=0; dy<dery+1; dy++){
        for (int dx=0; dx<derx+1; dx++){
            cout << derTensor({dx,dy}) << " " ;

        }cout << endl;
    }

    // =======================================================
    // Create pseudo mesh for testing  
    int M = 4;
    int N = 4;
    double hx = 1.0/(double)M;
    double hy = 1.0/(double)N;

    Tensor<vertex> mesh = Tensor<vertex>(2);

    mesh.setSize({M,N});

    for (int j=0; j<N; j++){
        for (int i=0; i<M; i++){
            mesh({i,j}) = {i*hx, j*hy};
        }
    }

    vector<vertex> cornerSet;
    cout << endl;
    cornerSet.push_back(mesh({1,1}));
    cout << mesh({1,1})[0] << "  " << mesh({1,1})[1] << endl;
    cornerSet.push_back(mesh({2,1}));
    cout << mesh({2,1})[0] << "  " << mesh({2,1})[1] << endl;
    cornerSet.push_back(mesh({2,2}));
    cout << mesh({2,2})[0] << "  " << mesh({2,2})[1] << endl;
    cornerSet.push_back(mesh({1,2}));
    cout << mesh({1,2})[0] << "  " << mesh({1,2})[1] << endl;
 
    degreex = 2;
    degreey = 2;

    polynomial testint = polynomial(degreex,degreey);

    double intcoef[4] = {0,0,0,1};

    testint.setCoef(intcoef);

    cout << "poly integral : "  << polyNumIntegralFace(cornerSet, 1.0, {0.0,0.0}, testint) << endl; 
    cout << "Previously defined integral : " << NumIntegralFace(cornerSet, {1,1}, {0.0,0.0}, 1.0, basePoly) << endl;

    // poly num matched with previously defined function

    return 0;
}
