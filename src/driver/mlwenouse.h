#ifndef MLWENOUSE_H_
#define MLWENOUSE_H_

#include "reconstruction.h"
#include "util.h"

/**!
 * An interface created to simplify usage of reconstruction.
 * Using MLWENO reconstruction should always be calling MLWENOUse class.
 * Instead of calling multiLevelReconstruction directly.
 */

namespace MLWENO{

    class MLWENOUse {
        public:
             //! A constructor
             /**! 
              * Construct a MLWENOUse class. 
              * This MLWENO containing multiple multi-level weno reconstruction instances.
              * Construct a MLWENOUse class by assigning a meshinfo pointer to it.
              * By using this MLWENOUse class, there is no need to separate boundary.
              * Different treatment can be applied to any place in the computational domain.
              */

             MLWENOUse() {};

             /**! A destructor
              * Destruct a MLWENOUse class.
              */
             
             ~MLWENOUse() {};
 
             /**
              * Create MLWENO reconstruction instance from a given MLWENOPrepare class.
              * We prepare for weno reconstruction only once.
              * We can have different MLWENO class, but only one WENOPrepare.
              * Create an MLWENO instance with given levels.
              */
             int AddMLWENOLevel(const std::string& location,
                                const unordered_set<std::string>& selectLevels,
                                MLWENOPrepare * mlpPtr);

             /**
              * Modify reconstruction stencils.
              * Assign ways to look for stencils.
              */
             void AssignWENOStencils(const std::string& location,
                                     std::string level, 
                                     vector<indice> newReconstMethod); 

             void AssignWENOStencils(const int& location,
                                     std::string level, 
                                     vector<indice> newReconstMethod); 

             /**
              * Evaluate a reconstruction value using defined MLWENO instances.
              */
             double Evaluate(const vertex& point, 
                             const indice& globalCell,
                             const MeshInfo& mi);
        private:

             /**
              * User-defined function.
              * Use this function to treat boundary differently.
              */
             int AssignInstance_(const indice& globalCell) const;

             /**
              * Location map.
              * Storing the information indicating how to assign MLWENO instance.
              */
             unordered_map<std::string, int> AssignMap_;

             vector<multiLevelReconstruction *> mlrIns_;
    };

}

#endif
