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

        double lambdad(const int& i,
                       const vertex& point)const;

        //! Rational function +- 1 on opposite edges
        //! Arbitrary values on other edges.
        double R(const int& e1,
                 const int& e2,
                 const vertex& point) const;

        //! Rational function with 1 on edge i
        //! 0 on the opposite edge
        //! Arbitrary values on other edges.
        double R(const int& e,
                 const vertex& point) const;

        //! Derivative of the rational auxilliary function
        vertex dR(const int& e1,
                  const int& e2,
                  const vertex& point) const;

        vertex dR(const int& e,
                  const vertex& point) const;

        double rational(const vertex& point) const;

        vertex dRational(const vertex& point) const;

        vertexSet corners() const {return corners_;}

        vertex unitnormal(const int& e) const {return unitNormals_.at(e);};

        //! Define lagrange basis polynomials
        std::array<double, 3> lagrangeE(const vertex& point, int nEdge, int jNode) const;
        std::array<double, 3> lagrangeV(const vertex& point, int nEdge, int i) const;

        double PhiSupp(const int& i, 
                       const vertex& point,
                       const int& r) const; 

    private:
        friend class BRMixed; 
        friend class Hdivmixed;

        int polynomial_degree_;

        //! Four corners of the given element
        vertexSet corners_;

        //! Four unit normal vectors of the corresponding edges
        vertexSet unitNormals_;

        //! Two unit normal vectors of the diagonal edges
        vertexSet unitNormals_d_;

        //! Four unit tangent vectors of the corresponding edges
        vertexSet unitTangents_;

        //! Storing edge nodes for lagrangian interpolation
        vertexSet lagEdgeNode_;

        //! Projection functions
        double projToEdge_(int Edge, const vertex& point) const;

        //! Calculate the distance between point and edge
        //! opposite to normal direction.
        double distance_(const vertexSet& edge, 
                         const vertex& point) const;
};

#endif
