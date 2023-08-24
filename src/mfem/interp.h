#ifndef INTERP_H_
#define INTERP_H_

#include "util.h"

namespace interpolation{
    class Lagrange {
        public:
            Lagrange(){};
            Lagrange(const tensor<int>& order);
            ~Lagrange(){clear();};

            /*!
             * Initiate a lagrange polynomial.
             * With given interpolation points.
             */
            void init(const tensorSet<double>& points);

            /*!
             * Initiate a lagrange polynomial. 
             * Without given interpolation points.
             * Default chebyshev interpolation points.
             */
            void init(const int order,
                      const tensorSet<double>& corners);

            /*!
             * Clear the defined lagrange polynomial.
             */
            void clear();

            /*!
             * Evaluate the defined lagrange polynomial.
             * Given a point.
             */
            double eval(const tensor<double>& point) const;

            /*!
             * Evaluate the derivative of the defined lagrange
             * polynomial. A point is given.
             */
            double evalDeriv(const tensor<double>& point) const;

        private:

            /*!
             * Polynomial order.
             */
            tensor<int> order_;

            /*!
             * Create chebyshev inerpolation points
             */
            vertexSet chebyshevPoints_(const vertex& left,  
                                       const vertex& right);

            /*!
             * Store interpolation points.
             */
            tensorSet<double> interpoPoints_;
    };
}

#endif
