#ifndef REACTION_H_
#define REACTION_H_

#include "reconstruction.h"

class reaction {
    public:
        reaction() {mlrPtr_ = new MLWENO::multiLevelReconstruction();};
        ~reaction() {delete mlrPtr_;};

    private:
        MLWENO::multiLevelReconstruction * mlrPtr_;
};

#endif
