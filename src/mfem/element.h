#ifndef ELEMENT_H_
#define ELEMENT_H_

#include "util.h"

/**
 * The class containing all the information
 * regarding element.
 * 3 - 2
 * |   |
 * 0 - 1
 * Vertex ordering as above
 */
class element{

    public:
        element();
        element(const MeshInfo& mi, const indice& global);
        ~element();

        //! Get four corners for this element
        void GetCorners(const vertexSet& corners);
        void GetCorners(const MeshInfo& mi, 
                        const indice& global){
            GetCorners(extractCorners(mi,global));
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
        double lambda(const int& e1, const int& e2,
                      const vertex& point) const;

    private:

        vertexSet corners_;

        // Calculate the distance between point and edge
        // opposite to normal direction.
        double distance_(const vertexSet& edge, 
                         const vertex& point) const;
};

#endif
