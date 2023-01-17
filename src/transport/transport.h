#ifndef TRANSPORT_H_
#define TRANSPORT_H_

#include <petsc.h>
#include <map>

#include "reconstruction.h"
#include "input.h"

namespace Transport{

    class advection {
        public:
            advection();
            ~advection() {Clear();};

            void Clear() const;

            double AdvFlux();

        private:
            std::map<indice, int> boundaryType_;
            std::map<indice, MLWENO::reconstruction *> advRecon_;

    };

    class diffusion {
        public:
            diffusion();
            ~diffusion() {Clear();};

            void Clear() const;

            double DiffFlux();

        private:
            std::map<indice, int> boundaryType_; 
            std::map<indice, MLWENO::reconstruction *> diffRecon_;

    };

    class transport : public advection, public diffusion{
        public:
            transport();
            ~transport() {};

            double TotalFlux();

        private:

    };

}
#endif
