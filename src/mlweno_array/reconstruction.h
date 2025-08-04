#ifndef RECONSTRUCTION_H_
#define RECONSTRUCTION_H_

#include "tensorstencilpoly.h"

// seriel version
class reconstruction {

    public:

        reconstruction()  {};
        ~reconstruction() {};

        int setWgts(vector<> ){

            return 1;
        }

    private:

        vector<double> linwgts;
        vector<double> nonlinwgts;

        vector<tensorstencilpoly *> stenpoly;

}

#endif
