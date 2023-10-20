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

        //! Constant part
        double phic(const basis& auxilliary,
                    const int& nEdge, 
                    const vertex& point) const;

        //! Linear part
        double phil(const basis& auxilliary,
                    const int& nEdge,
                    const vertex& point) const;

        // ! Test function of H(div) mixed function space
        void Test();
   

    private:
        

};

#endif
