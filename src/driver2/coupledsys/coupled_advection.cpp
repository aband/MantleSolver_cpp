#include "driver.h"

int Driver::computephase(const std::vector<vertex>& gaussp,
                         const vertexSet& edgep,
                         const indice& gcell,
                         const Tensor<weights>& allwgtsHD, double ** lHD,
                         const Tensor<weights>& allwgtsCD, double ** lCD,
                         vector<double>& kappa,
                         vector<double>& lambda,
                         vector<double>& phif){

    // Compute phase attributes on both sides of the edge
    // Perform harmonic average of two different values
    double HD = 0.0, CD = 0.0;
    for (int g=0; g<gaussp.size(); g++){
        HD = advection.eval(gaussp.at(g), ml, "all", 
             allwgtsHD({gcell[0], gcell[1]}), gcell, lHD);

        CD = advection.eval(gaussp.at(g), ml, "all", 
             allwgtsCD({gcell[0], gcell[1]}), gcell, lCD);

        // Calculate volumetric fraction at given quadrature points
        double depth  = myPhase->pPtr->GetDepth(gaussp.at(g)[1], myPhase->pp->l0);  
        double lithoP = myPhase->pPtr->GetScaledLithoP(depth);

        myPhase->pPtr->evalPhase(HD,CD,lithoP);
        phif.at(g) = myPhase->pPtr->phi.mlt;

        kappa.at(g) = phif.at(g) * myPhase->pPtr->Getcf(CD) / CD;
        lambda.at(g) = phif.at(g) * myPhase->pPtr->Getef() / HD; 
    }

    return 1;
}

inline vertex effectVel(const vertex& vf, const vertex& vs, const double& c, const double& phif){

    return vf*phif*c + (1-c)*vs;
}

// Compute effective velocity (interior edges)
// Two cell indice should be provided
int Driver::computeEffVel(const vector<vertex>& gaussp,
                          const vertexSet& edgep,
                          const indice& gcellin, const indice& gcellout,
                          const Tensor<weights>& allwgtsHD, double ** lHD,
                          const Tensor<weights>& allwgtsCD, double ** lCD,
                          vector<vertex>& effvelHD, 
                          vector<vertex>& effvelCD){

    vector<double> phifin ; phifin.resize(gaussp.size());
    vector<double> phifout; phifout.resize(gaussp.size());
    vector<double> lambdain ; lambdain.resize(gaussp.size());
    vector<double> lambdaout; lambdaout.resize(gaussp.size());
    vector<double> kappain ; kappain.resize(gaussp.size());
    vector<double> kappaout; kappaout.resize(gaussp.size());

    // Extract velocity on these given gauss points
    vector<vertex> vel_darcy = 
    ExtractVelocity(&sresult_->vel_darcy, &sresult_->g_darcy,
                    refArrayDarcyEssen_,mi,
                    gaussp, gcellin,*hdiv_,*basis_,{1});
    
    vector<vertex> vel_stokes = 
    ExtractVelocity(&sresult_->vel_stokes, &sresult_->g_stokes,
                    refArrayStokesEssen_,mi,
                    gaussp, gcellin,*br_,*basis_,{1});

    computephase(gaussp, edgep, gcellin, allwgtsHD, lHD, allwgtsCD, lCD,
                 kappain, lambdain, phifin);

    computephase(gaussp, edgep, gcellout, allwgtsHD, lHD, allwgtsCD, lCD,
                 kappaout, lambdaout, phifout);

    for (int g=0; g<gaussp.size(); g++){
        double kappa_mean = harmonic_mean(kappain.at(g), kappaout.at(g));
        double lambda_mean = harmonic_mean(lambdain.at(g), lambdaout.at(g));

        double phif_mean = harmonic_mean(phifin.at(g), phifout.at(g));
        effvelHD.at(g) = effectVel(vel_darcy.at(g), vel_stokes.at(g), lambda_mean, phif_mean);
        effvelCD.at(g) = effectVel(vel_darcy.at(g), vel_stokes.at(g), kappa_mean , phif_mean);
    }

    return 1;
} 

// Used on the boundary where only one sided value will be used
int Driver::computeEffVel(const vector<vertex>& gaussp,
                          const vertexSet& edgep,
                          const indice& gcell,
                          const Tensor<weights>& allwgtsHD, double ** lHD,
                          const Tensor<weights>& allwgtsCD, double ** lCD,
                          vector<vertex>& effvelHD, 
                          vector<vertex>& effvelCD){

    vector<double> phif ; phif.resize(gaussp.size());
    vector<double> lambda ; lambda.resize(gaussp.size());
    vector<double> kappa ; kappa.resize(gaussp.size());

    // Extract velocity on these given gauss points
    vector<vertex> vel_darcy = 
    ExtractVelocity(&sresult_->vel_darcy, &sresult_->g_darcy,
                    refArrayDarcyEssen_,mi,
                    gaussp, gcell,*hdiv_,*basis_,{1});
    
    vector<vertex> vel_stokes = 
    ExtractVelocity(&sresult_->vel_stokes, &sresult_->g_stokes,
                    refArrayStokesEssen_,mi,
                    gaussp, gcell,*br_,*basis_,{1});

    computephase(gaussp, edgep, gcell, allwgtsHD, lHD, allwgtsCD, lCD,
                 kappa, lambda, phif);

    for (int g=0; g<gaussp.size(); g++){

        effvelHD.at(g) = effectVel(vel_darcy.at(g), vel_stokes.at(g), lambda.at(g), phif.at(g));
        effvelCD.at(g) = effectVel(vel_darcy.at(g), vel_stokes.at(g), kappa.at(g) , phif.at(g));
    }

    return 1;
}

int Driver::updateEdgeFlux(Tensor<double>& vertedgeHD, Tensor<double>& horiedgeHD,
                           Tensor<double>& vertedgeCD, Tensor<double>& horiedgeCD,
                           const Tensor<weights>& allwgtsHD, double ** lHD,
                           const Tensor<weights>& allwgtsCD, double ** lCD){

    Tensor_zero(vertedgeHD);
    Tensor_zero(horiedgeHD);
    Tensor_zero(vertedgeCD);
    Tensor_zero(horiedgeCD);

    const valarray<double>& gwe = GaussWeightsEdge;
    const valarray<double>& gpe = GaussPointsEdge;

    std::vector<vertex> gaussp;
    gaussp.resize(gpe.size());

    for (int j=0; j<mi.MPIglobalCellSize[1]; j++){
    for (int i=0; i<mi.MPIglobalCellSize[0]; i++){

        double flux = 0.0;
        indice globalcell {i,j};
        indice cellout;

        // Extract corners with respect to given global indice
        vertexSet corners = extractCorners(mi, globalcell); 

        // =============================================================
        // Get horizontal edge
        vertexSet hori {corners.at(0), corners.at(1)};
      
        // Extract velocity on this edge 
        for (int g=0; g<gpe.size(); g++){
            gaussp.at(g) = GaussMapPointsEdge({gpe[g]},hori);
        }   

        if (j==0){
           //bottom edge, fix values

        } else {

            cellout = globalcell + mi.faceNormal[3];
        }

        horiedgeHD({i,j}) = ;
        horiedgeCD({i,j}) = ;

        // =============================================================
        // Get vertical edge
        vertexSet vert {corners.at(3), corners.at(0)};

        gaussp.clear();gaussp.resize(gpe.size());
        // Extract velocity on this edge 
        for (int g=0; g<gpe.size(); g++){
            gaussp.at(g) = GaussMapPointsEdge({gpe[g]},hori);
        }   

        // boundary
        if (i==0){
            // no flow boundary on left side
            flux    = 0.0;
        } else {
            cellout = globalcell + mi.faceNormal[3];

        }
        vertedgeHD({i,j}) = ;
        vertedgeHD({i,j}) = ;

    }}

    // right side i = mi.MPIglobalCellSize[0];
    // no flow flux = 0;


    return 1;
}
