#ifndef BRMIXED_H_
#define BRMIXED_H_

#include "basis.h" 

// An enriched DS1 H1 conforming finite element space

class BRMixed {

    public:

        BRMixed() {name = "BR"; elemDOF = 12;};
        ~BRMixed() {};

        //! Mapping from local degree of freedom to
        //! global index of degree of freedom
        std::array<int,12> LocalToGlobal(const MeshInfo& mi,
                                         const indice& globalElement) const; 

        void ComputeTotalDOF(const MeshInfo& mi);

        std::string Name() const{return name;}

        int getDOF() const {return totalDOF_;};

        //! Nodal basis function
        double phiv(const basis& basis_,
                    const int& nnodal,
                    const vertex& point) const;

        vertex dPhiv(const basis& basis_,
                     const int& nnodal,
                     const vertex& point) const;

        //! Bubble function
        double phie(const basis& basis_,
                    const int& nEdge,
                    const vertex& point) const;

        vertex dPhie(const basis& basis_,
                     const int& nEdge,
                     const vertex& point)const;

        void Test(const basis& basis_);

        std::array<std::array<double,4>,12> ComputeGradBRmixed(const basis& basis_,
                                                               const vertex& point) const;
        // ! Return all basis functions evaluated at a given point
        std::array<vertex, 12> ComputeBRmixed(const basis& basis_,
                                              const vertex& point) const;

        // ! Return the basis function on the given edge evaluated at a given point
        vertex ComputeBRmixed(const basis& basis_, 
                              const vertex& point,
                              const int& local) const;

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

        // ! Return global cell index, corresponding local index in this order
        std::vector<int> GlobalToLocalMapBndry(const MeshInfo& mi,
                                               const int& gdof) const;

        bool onBndry(const MeshInfo& mi,
                     const int& globaldof)const;

        double Pressure() const {return 1.0;};

        std::string name;

        int elemDOF;

        double R(const basis& basis_,
                  const vertex& point) const{
        return R_(basis_, point);
        }


    private:

        double phie_(const basis& basis_,
                     const int& nEdge,
                     const vertex& point) const;

        double R_(const basis& basis_,
                  const vertex& point) const;

        vertex dR_(const basis& basis_,
                   const vertex& point) const;

        vertex dphie_(const basis& basis_,
                      const int& nEdge,
                      const vertex& point) const;

        int totalDOF_;

};

#endif
