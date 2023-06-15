#ifndef HDIVMIXED_H_
#define HDIVMIXED_H_

#include "basis.h"

class Hdivmixed: virtual public basis{
    public: 
        Hdivmixed() {};
        ~Hdivmixed() {};

        // ! Test function of H(div) mixed function space
        void Test();
    
};

#endif
