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

    vector<vertex> effvelHD; effvelHD.resize(gaussp.size());
    vector<vertex> effvelCD; effvelCD.resize(gaussp.size());

    for (int j=0; j<mi.MPIglobalCellSize[1]; j++){
    for (int i=0; i<mi.MPIglobalCellSize[0]; i++){

        double fluxHD = 0.0;
        double fluxCD = 0.0;
        indice gcell {i,j};
        indice cellout;
        effvelHD.clear(); effvelHD.resize(gaussp.size());
        effvelCD.clear(); effvelCD.resize(gaussp.size());

        // Extract corners with respect to given global indice
        vertexSet corners = extractCorners(mi, gcell); 

        // =============================================================
        // Get horizontal edge
        vertexSet hori {corners.at(0), corners.at(1)};
      
        // Extract velocity on this edge 
        for (int g=0; g<gpe.size(); g++){
            gaussp.at(g) = GaussMapPointsEdge({gpe[g]},hori);
        }   

        if (j==0){
            // Use one sided velocity
            computeEffVel(gaussp, hori, gcell, allwgtsHD, lHD, allwgtsCD, lCD, effvelHD, effvelCD);
            // bottom edge, fix values
            fluxHD = edgefluxintegral(hori, HDbottom, effvelHD);
            fluxCD = edgefluxintegral(hori, CDbottom, effvelCD);
        } else {

            cellout = gcell + mi.faceNormal[0];
            computeEffVel(gaussp, hori, gcell, cellout, allwgtsHD, lHD, allwgtsCD, lCD, effvelHD, effvelCD);
 
            fluxHD = edgefluxintegral(mi, gcell, cellout, hori, allwgtsHD, effvelHD, ml, advection, lHD);
            fluxCD = edgefluxintegral(mi, gcell, cellout, hori, allwgtsCD, effvelCD, ml, advection, lCD);
        }

        horiedgeHD({i,j}) = fluxHD;
        horiedgeCD({i,j}) = fluxCD;

        // =============================================================
        // Get vertical edge
        vertexSet vert {corners.at(3), corners.at(0)};

        gaussp.clear();gaussp.resize(gpe.size());
        // Extract velocity on this edge 
        for (int g=0; g<gpe.size(); g++){
            gaussp.at(g) = GaussMapPointsEdge({gpe[g]},hori);
        }   

        effvelHD.clear(); effvelHD.resize(gaussp.size());
        effvelCD.clear(); effvelCD.resize(gaussp.size());

        // boundary
        if (i==0){
            // no flow boundary on left side
            fluxHD    = 0.0;
            fluxCD    = 0.0;
        } else {
            cellout = gcell + mi.faceNormal[3];
            computeEffVel(gaussp, vert, gcell, cellout, allwgtsHD, lHD, allwgtsCD, lCD, effvelHD, effvelCD);

            fluxHD = edgefluxintegral(mi, gcell, cellout, vert, allwgtsHD, effvelHD, ml, advection, lHD);
            fluxCD = edgefluxintegral(mi, gcell, cellout, vert, allwgtsCD, effvelCD, ml, advection, lCD);
        }
        vertedgeHD({i,j}) = fluxHD;
        vertedgeCD({i,j}) = fluxCD;

    }}

    // right side i = mi.MPIglobalCellSize[0];
    // no flow flux = 0;
    // no additional operations needed
    
    // On top side j = mi.MPIglobalCellSize[1]
    // Free outflow
    for (int i=0; i<mi.MPIglobalCellSize[0]; i++){

        // Regarded as outside cell
        indice gcell {i, mi.MPIglobalCellSize[1]-1};
        vertexSet corners = extractCorners(mi, gcell);
        vertexSet hori    = {corners.at(3), corners.at(2)};

        for (int g=0; g<gpe.size(); g++){gaussp.at(g) = GaussMapPointsEdge({gpe[g]}, hori);}

        computeEffVel(gaussp, hori, gcell, allwgtsHD, lHD, allwgtsCD, lCD, effvelHD, effvelCD);

        horiedgeHD({i, mi.MPIglobalCellSize[1]}) = edgefluxintegral(mi, gcell, hori, allwgtsHD, effvelHD, ml, advection, lHD);
        horiedgeCD({i, mi.MPIglobalCellSize[1]}) = edgefluxintegral(mi, gcell, hori, allwgtsCD, effvelCD, ml, advection, lCD);
    }

    return 1;
}

int Driver::getflux(const Tensor<weights>& allwgtsHD, double ** lHD, 
                    const Tensor<weights>& allwgtsCD, double ** lCD, 
                    double **lfHD, double** lfCD){

    Tensor<double> horiedgefluxHD = Tensor<double>(2);
    horiedgefluxHD.setSize({mi.MPIlocalCellSize[0], mi.MPIlocalCellSize[1]+1});

    Tensor<double> vertedgefluxHD = Tensor<double>(2);
    vertedgefluxHD.setSize({mi.MPIlocalCellSize[0]+1, mi.MPIlocalCellSize[1]});

    Tensor<double> horiedgefluxCD = Tensor<double>(2);
    horiedgefluxCD.setSize({mi.MPIlocalCellSize[0], mi.MPIlocalCellSize[1]+1});

    Tensor<double> vertedgefluxCD = Tensor<double>(2);
    vertedgefluxCD.setSize({mi.MPIlocalCellSize[0]+1, mi.MPIlocalCellSize[1]});

    updateEdgeFlux(vertedgefluxHD, horiedgefluxHD,
                   vertedgefluxCD, horiedgefluxCD,
                   allwgtsHD, lHD,
                   allwgtsCD, lCD);

    for (int j=0; j<mi.MPIglobalCellSize[1]; j++){
    for (int i=0; i<mi.MPIglobalCellSize[0]; i++){

        lfHD[j][i] = getcellflux(mi, {i,j}, vertedgefluxHD, horiedgefluxHD);
        lfCD[j][i] = getcellflux(mi, {i,j}, vertedgefluxCD, horiedgefluxCD);

    }}

    return 1;
}
