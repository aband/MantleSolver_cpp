#ifndef MLWENO_H_
#define MLWENO_H_

#include "reconstruction.h"
#include "util.h"

/**!
 * An interface created to simplify usage of reconstruction.
 * Using MLWENO should always be calling this class.
 * Instead of calling multiLevelReconstruction directly.
 */

namespace MLWENO{

    class MLWENO {
        public:
             //! A constructor
             /**! 
              * Construct a MLWENO class. 
              * This MLWENO containing multiple multi-level weno reconstruction instances.
              */

             MLWENO() {};

             /**! A destructor
              * Destruct a MLWENO class.
              */
             
             ~MLWENO() {};
 
             /**
              * Create MLWENO reconstruction instance from a given MLWENOPrepare class.
              * We prepare for weno reconstruction only once.
              * We can have different MLWENO class, but only one WENOPrepare.
              * Create an MLWENO instance with given levels.
              */

             void AddMLWENOInstance(const vector<std::string>& selectLevels,
                                    MLWENOPrepare * mlpPtr);

        private:

             vector<multiLevelReconstruction * mlrPtr> mlrIns;


    };

}

#endif
