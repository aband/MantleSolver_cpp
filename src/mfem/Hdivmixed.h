#ifndef HDIVMIXED_H_
#define HDIVMIXED_H_

#include "basis.h"

// First order H(div) conforming mixed space.
// Derived from Direct Serendipity space. a BDM style element.
// Basis functions are constructed based on method
// mentioned in Direct Serendipity space.
class Hdivmixed{
    public: 
        Hdivmixed() {};
        ~Hdivmixed() {};

        std::array<int, 8>  LocalToGlobal(const MeshInfo& mi,
                                          const indice& globalElement) const;

        //! Constant part
        vertex phic(const basis& basis_,
                    const int& nEdge, 
                    const vertex& point) const;

        //! Linear part
        vertex phil(const basis& basis_,
                    const int& nEdge,
                    const vertex& point) const;

        void ComputeTotalDOF(const MeshInfo& mi);

        int getDOF() const {return totalDOF_;};

        //! Divergence of the constant part
        double divphic(const basis& basis_,
                       const int& nEdge,
                       const vertex& point) const; 

        // ! Test function of H(div) mixed function space
        void Test(const basis& basis_);
 
        // ! Return all basis function evaluated at a given point
        std::array<vertex,8> ComputeHdivmixed(const basis& basis_,
                                              const vertex& point) const;

        // ! Return the basis function on the given edge evaluated at a given point
        std::array<vertex,2> ComputeHdivmixed(const basis& basis_,
                                              const vertex& point,
                                              const int& edge) const;

    private:
       
        vertex curlLambda_(const basis& basis_,
                           const int& e) const; 

        vertex curlR_(const basis& basis_,
                      const int& e1,
                      const int& e2,
                      const vertex& point) const;

        vertex curlR_(const basis& basis_,
                      const int& e,
                      const vertex& point) const;

        double phie_(const basis& basis_,
                     const int& e,
                     const vertex& point) const;

        double phiv_(const basis& basis_,
                     const int& e,
                     const vertex& point) const;

        vertex curlphiv_(const basis& basis_,
                         const int& e,
                         const vertex& point) const;

        vertex curlphiv_star_(const basis& basis_,
                              const int& e,
                              const vertex& point) const;

        vertex phiv_star_star_(const basis& basis_,
                               const int& i,
                               const vertex& point) const;

        // Not normalized yet
        vertex phil_(const basis& basis_,
                     const int& nEdge,
                     const vertex& point) const;

        int totalDOF_;
};

#endif
