#ifndef MLWENOUSE_H_
#define MLWENOUSE_H_

#include "reconstMLWENO.h"
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
 
             /**!
              * Create MLWENO reconstruction instance from a given MLWENOPrepare class.
              * We prepare for weno reconstruction only once.
              * We can have different MLWENO class, but only one WENOPrepare.
              * Create an MLWENO instance with given levels.
              */
             int AddMLWENOLevel(const std::string& location,
                                const unordered_set<std::string>& selectLevels,
                                MLWENOPrepare * mlpPtr);

             /**!
              * Modify reconstruction stencils.
              * Assign ways to look for stencils.
              */
             void AssignWENOStencils(const std::string& location,
                                     const std::string& level, 
                                     const vector<indice>& newReconstMethod); 

             void AssignWENOStencils(const int& location,
                                     const std::string& level, 
                                     const vector<indice>& newReconstMethod); 

             /**!
              * Assing linear weights to corresponding reconstruction levels. 
              * No need to be able to sum up to 1.
              */
             void AssignLinearWgts(const std::string& location,
                                   const std::string& level,
                                   const vector<double>& linWgts);

             void AssignLinearWgts(const int& location,
                                   const std::string& level,
                                   const vector<double>& linWgts);

             /**!
              * Update non linear weights.
              */
             void UpdateNonLinearWgts(const MeshInfo& mi, 
                                      const std::string& location,
                                      const std::string& weightType,
                                      bool (*func)(const indice& globalCell,
                                                   const MeshInfo& mi));

             void UpdateNonLinearWgts(const MeshInfo& mi, 
                                      const int& location,
                                      const std::string& weightType,
                                      bool (*func)(const indice& globalCell,
                                                   const MeshInfo& mi));

             /**!
              * Updated version.
              * No need of weighting type
              */
             void UpdateNonLinearWgts(const MeshInfo& mi, 
                                      const std::string& location,
                                      bool (*func)(const indice& globalCell,
                                                   const MeshInfo& mi));

             void UpdateNonLinearWgts(const MeshInfo& mi, 
                                      const int& location,
                                      bool (*func)(const indice& globalCell,
                                                   const MeshInfo& mi));

             void UpdateNonLinearWgts(const MeshInfo& mi, 
                                      const std::string& location,
                                      bool (*func)(const indice& globalCell,
                                                   const MeshInfo& mi),
                                      const std::string& name);

             void UpdateNonLinearWgts(const MeshInfo& mi, 
                                      const int& location,
                                      bool (*func)(const indice& globalCell,
                                                   const MeshInfo& mi),
                                      const std::string& name);

             /**!
              * Evaluate a reconstruction value using defined MLWENO instances.
              */
             double Evaluate(const vertex& point, 
                             const indice& globalCell,
                             const MeshInfo& mi,
                             const int& location) const;

             double Evaluate(const vertex& point,
                             const indice& globalCell,
                             const MeshInfo& mi,
                             const std::string& location) const;

             /**!
              * Print non-linear weights relating to weights.
              */
             void PrintNonLinearWgts(const int& location,
                                    const MeshInfo& mi);

             void PrintNonLinearWgts(const std::string& location,
                                    const MeshInfo& mi);

        private:

             /**
              * Location map.
              * Storing the information indicating how to assign MLWENO instance.
              */
             unordered_map<std::string, int> AssignMap_;

             /**!
              * Vector holding pointers to multi level reconstructions
              */
             vector<multiLevelReconstruction *> mlrIns_;
    };

}

#endif
