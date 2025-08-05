#ifndef RECONSTRUCTION_H_
#define RECONSTRUCTION_H_

#include "tensorstencilpoly.h"

// seriel version
class reconstruction {

    public:

        reconstruction()  {};
        ~reconstruction() {};

        int setWgts(vector<vector<sigma>>& );

    private:

        vector<int> stencilnum;

        vector<vector<double>> linwgts;
        vector<vector<double>> nonlinwgts;

        vector<vector<indice>> pickstencil;
}

const int geteta(int r){

    if (r==0){
        return 1;
    } else if (r==1){
        return 3;
    } else {
        return 4;
    }
}

#endif
