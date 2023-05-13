#ifndef DIFFUSION_H_
#define DIFFUSION_H_

#include "reconstruction.h"

namespace SymDiffusion{
    class diffusion {
        public:
            diffusion() {mlrPtrHori_ = new MLWENO::multiLevelReconstruction();
                         mlrPtrVert_ = new MLWENO::multiLevelReconstruction();};
            ~diffusion() {delete mlrPtrHori_; delete mlrPtrVert_;};

            /**
             * Diffusion defined on edges.
             * Different reconstruction methods used for horizontal and veritical edges.
             */
            unordered_map<std::string, vector<indice>> reconstMethodsVert;
            unordered_set<std::string> wenoLevelsVert;

            unordered_map<std::string, vector<indice>> reconstMethodsHori;
            unordered_set<std::string> wenoLevelsHori;

            // Interpolation positions
            const double alpha = 0.5;
            const double beta = 1.5;

            // Create corresponding multilevel reconstruction with the given information.
            void CreateMLWENO(const MLWENO::MLWENOPrepare& mlp, const MeshInfo& mi);

            // Update nonlinear weights
            void UpdateNonLinearWgts(const MeshInfo& mi);

            // Update flux and its derivative across the edges collectively
            void UpdateEdgeFlux(const MeshInfo& mi);

            void UpdateEdgeFluxDerivative(const MeshInfo& mi);

            // Return sum of integrated flux or corresponding derivative for a given target cell
            double Flux(const MeshInfo& mi, const indice& global);

            unordered_map<int, double> derivFlux(const MeshInfo& mi, const indice& global);

            // Print information
            void GetInfo(const MeshInfo& mi);

        private:

            double fixed_ = 0.0;

            // Calculate diffusive flux.
            // And its corresponding derivatives.
            double flux_(const double * ru, int n);

            double dflux_(const double * ru, int n);

            //! Define boundary conditions here
            double boundaryCondition_(double * ru, int n, const int& flag);

            //! Compute integrated diffusion flux on the given edge
            //! Horizontal or vertical multilevel reconstruction is passed in
            //! as a parameter.
            double edgeFlux_(const MeshInfo& mi, 
                             const indice& globalCell,
                             const int& k,
                             const int& size,
                             const vertexSet& edge,
                             const double& scale,
                             const MLWENO::multiLevelReconstruction& mlrPtr);

            unordered_map<int, double> derivEdgeFlux_(const MeshInfo& mi,
                                                      const indice& globalCell,
                                                      const int& k,
                                                      const int& size,
                                                      const vertexSet& edge,                       
                                                      const double& scale,
                                                      const MLWENO::multiLevelReconstruction& mlrPtr);

            /**
             * Determine whether the target cell is inside the boundary layer
             */
            int Interior_(const int& k, const int& size);
    
            //! Store calculated diffusive flux on the edge
            //! And its corresponding derivative
            unordered_map<int, double> edgeHoriFlux_;
            unordered_map<int, double> edgeVertFlux_;

            unordered_map<int, unordered_map<int, double>> derivEdgeHoriFlux_;
            unordered_map<int, unordered_map<int, double>> derivEdgeVertFlux_;

            /**
             * Pointer to a multilevel reconstruction class.
             */
            MLWENO::multiLevelReconstruction * mlrPtrHori_;
            MLWENO::multiLevelReconstruction * mlrPtrVert_;
};
}

namespace NonSymDiffusion{
    class diffusion {
        public:
            diffusion() {mlrPtrHoriUp_    = new MLWENO::multiLevelReconstruction();
                         mlrPtrHoriDown_  = new MLWENO::multiLevelReconstruction();
                         mlrPtrVertRight_ = new MLWENO::multiLevelReconstruction();
                         mlrPtrVertLeft_  = new MLWENO::multiLevelReconstruction();};

            ~diffusion() {delete mlrPtrHoriUp_; 
                          delete mlrPtrHoriDown_; 
                          delete mlrPtrVertRight_; 
                          delete mlrPtrVertLeft_;};

            unordered_map<std::string, vector<indice>> reconstMethodsVertRight;
            unordered_map<std::string, vector<indice>> reconstMethodsVertLeft;
            unordered_map<std::string, vector<indice>> reconstMethodsHoriUp;
            unordered_map<std::string, vector<indice>> reconstMethodsHoriDown;

            unordered_set<std::string> wenoLevelsVertRight;
            unordered_set<std::string> wenoLevelsVertLeft;
            unordered_set<std::string> wenoLevelsHoriUp;
            unordered_set<std::string> wenoLevelsHoriDown;

            // Interpolation positions
            const double alpha = 0.5;
            const double beta = 1.5;

            void CreateMLWENO(const MLWENO::MLWENOPrepare& mlp, const MeshInfo& mi);

            void UpdateNonLinearWgts(const MeshInfo& mi);

            void UpdateEdgeFlux(const MeshInfo& mi);

            double Flux(const MeshInfo& mi, const indice& global);

            void GetInfo(const MeshInfo& mi);

        private:
            double flux_(const std::array<double,4>& ru);
            double dflux_(const std::array<double,4>& dru);

            double boundaryCondition_(std::array<double,4>& ru, const std::array<int,2>& posOut);

            bool Interior_(const MeshInfo& mi, const indice& target);

            double edgeFlux_(const MeshInfo& mi,
                             const indice& globalIn,
                             const indice& globalOut,
                             const vertexSet& edge,
                             const double& scale,
                             const std::array<int, 2>& posIn,
                             const std::array<int, 2>& posOut,
                             const MLWENO::multiLevelReconstruction& mlrPtrIn,
                             const MLWENO::multiLevelReconstruction& mlrPtrOut);

            //! Store calculated diffusive flux on the edge
            //! And its corresponding derivative
            unordered_map<int, double> edgeHoriFlux_;
            unordered_map<int, double> edgeVertFlux_;

            /**
             * Pointers to different multilevel reconstruction class
             */
            MLWENO::multiLevelReconstruction * mlrPtrHoriUp_;
            MLWENO::multiLevelReconstruction * mlrPtrHoriDown_;
            MLWENO::multiLevelReconstruction * mlrPtrVertRight_;
            MLWENO::multiLevelReconstruction * mlrPtrVertLeft_;
    };
}

#endif
