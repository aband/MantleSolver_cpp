#ifndef HDIVMIXED_H_
#define HDIVMIXED_H_

#include "basis.h"

// First order H(div) conforming mixed space.
// Derived from Direct Serendipity space. a BDM style element.
// Basis functions are constructed based on method
// mentioned in Direct Serendipity space.
class Hdivmixed{
    public: 
        Hdivmixed() {name = "BDM";elemDOF = 8;};
        ~Hdivmixed() {};

        std::array<int, 8> LocalToGlobal(const MeshInfo& mi,
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

        std::string Name() const{return name;}

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

        // ! Standarized function usage
        std::vector<vertex> EvaluateAll(const basis& basis_,
                                        const vertex& point) const;

        vertex Evaluate(const basis& basis_,
                        const vertex& point,
                        const int& localdof) const;

        std::vector<std::array<double,4>> EvaluateGradAll(const basis& basis_,
                                                          const vertex& point) const;

        std::vector<int> LocalGlobalMap(const MeshInfo& mi,
                                        const indice& global) const;

        std::vector<int> GlobalToLocalMapBndry(const MeshInfo& mi,
                                               const int& globaldof) const;

        bool onBndry(const MeshInfo& mi,
                     const int& globaldof) const;

        double Pressure() const{return 1.0;};

        std::string name;

        int elemDOF;

        double phie(const basis& basis_,
                    const int& e,
                    const vertex& point) const{
        return phie_(basis_,e,point);
		  }
        double phiv(const basis& basis_,
                    const int& e,
                    const vertex& point) const{
        return phiv_(basis_,e,point);
		  }


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
