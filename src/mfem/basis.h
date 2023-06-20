#ifndef BASIS_H_
#define BASIS_H_

#include "util.h"

/*!
 * The class containing all the information
 * regarding element.
 * 3 - 2
 * |   |
 * 0 - 1
 * Vertex ordering as above.
 *
 * |-3-|
 * 0   2
 * |-1-|
 * Edge ordering as above.
 *
 * Class basis is used by assigning values
 * to the reference of four corners as a quadrilateral.
 */
class basis{

    public:
        basis() {};
        ~basis() {};

        //! Get four corners for this element
        void GetCorners(const vertexSet& corners);
        void GetCorners(const MeshInfo& mi, 
                        const indice& global){
            GetCorners(extractCorners(mi,global));
        }

        //! Overload ostream function
        friend std::ostream& operator<<(std::ostream& os, const basis& obj){
            os << "The element has four corners: " << endl;
            os << std::setw(3) << std::setprecision(3);
            for (const auto& it: obj.corners_){
                os << "( " << it[0] << ", " << it[2] << " )  ";
            }
            os << endl;
            return os;
        }

        // ===== Test =====
        void Test(const vertex& point);

    protected:

        //! Defined linear polynomial giving the distance
        //! between point and edge opposite the normal
        //! direction.
        double lambda(const int& e, 
                      const vertex& point) const;

        //! Defined linear polynomial giving the distance
        //! between point and diagonal opposite t
        //! e1 and e2 choosing 0,2 or 1,3
        double lambda(const int& e1, 
                      const int& e2,
                      const vertex& point) const;

        //! Rational function +- 1 on opposite edges
        //! Arbitrary values on other edges.
        double R(const int& e1,
                 const int& e2,
                 const vertex& point) const;

        //! Rational function with 1 on edge i
        //! 0 ont the opposite edge
        //! Arbitrary values on other edges.
        double R(const int& e,
                 const vertex& point) const;

        //! Define lagrange basis polynomials
        double lagrange(const vertex& point) const;

        //! Basis functions on cell, edges and vertex
        //! Satisfying unisolvence condition.
        bool Phi(const int i,
                 const int j) const;

        double phi() const;


    private:

        vertexSet corners_;

        //! Calculate the distance between point and edge
        //! opposite to normal direction.
        double distance_(const vertexSet& edge, 
                         const vertex& point) const;
};

#endif
