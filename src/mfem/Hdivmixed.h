#ifndef HDIVMIXED_H_
#define HDIVMIXED_H_

#include "basis.h"

// Second order H(div) conforming mixed space.
// Basis functions are constructed based on method
// mentioned in Direct Serendipity space.
class Hdivmixed{
    public: 
        Hdivmixed() {};
        ~Hdivmixed() {};


        // ! Test function of H(div) mixed function space
        void Test();
    
};

#endif
