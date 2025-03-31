#include "driver.h"

inline int clearall(vector<vertex>& effvel,
                    vector<vertex>& phasevel,
                    vector<vertex>& solidvel,
                    vector<double>& TDin,
                    vector<double>& TDout,
                    vector<double>& dTdHin,
                    vector<double>& dTdHout,
						  vector<double>& CDin,
						  vector<double>& CDout,
						  vector<double>& HDin,
						  vector<double>& HDout){

    int size = effvel.size();

    effvel.clear(); effvel.resize(size);
    phasevel.clear(); phasevel.resize(size);
    solidvel.clear(); solidvel.resize(size);
    TDin.clear(); TDin.resize(size);
    TDout.clear(); TDout.resize(size);
    dTdHin.clear(); dTdHin.resize(size);
    dTdHout.clear(); dTdHout.resize(size);
    CDin.clear(); CDin.resize(size);
    CDout.clear(); CDout.resize(size);
    HDin.clear(); HDin.resize(size);
    HDout.clear(); HDout.resize(size);

    return 1;
}

inline bool notoutflow(const double& flux){

    bool flag = false;

    if (flux < 0){
        flag = true;
    }

    return flag;
}

// Compute phase values on a given edge with respect to 
// a given cell
int Driver::computephase(const std::vector<vertex>& gaussp,
                         const vertexSet& edgep,
                         const indice& gcell,
                         const Tensor<weights>& allwgtsHD, double ** lHD,
                         const Tensor<weights>& allwgtsCD, double ** lCD,
                         vector<double>& cs,
                         vector<double>& cl,
                         vector<double>& phi,
                         vector<double>& TD,
                         vector<double>& dTdH,
	                      vector<double>& cd,
                         vector<double>& hd){

    double HD = 0.0, CD = 0.0;

    for (int g=0; g<gaussp.size(); g++){
        HD = advection.eval(gaussp.at(g), ml, location(mi, gcell), 
             allwgtsHD({gcell[0], gcell[1]}), gcell, lHD);

        CD = advection.eval(gaussp.at(g), ml, location(mi, gcell), 
             allwgtsCD({gcell[0], gcell[1]}), gcell, lCD);

        hd.at(g) = HD;
        cd.at(g) = CD;

        // Evaluate phase with given pressure, bulk enthalpy and bulk composition 
        double lithoP = myPhase->pPtr->GetStaticP(-1*gaussp.at(g)[1], 
                                                  myPhase->pPtr->l0); 
        myPhase->pPtr->evalPhase(HD, CD, lithoP);
        phi.at(g) = myPhase->pPtr->pc.phil;
        cs.at(g)  = myPhase->pPtr->pc.phi2/ (myPhase->pPtr->pc.phi2 +
                                             myPhase->pPtr->pc.phi1);
        cl.at(g)  = myPhase->pPtr->pc.cl;

        TD.at(g)  = myPhase->pPtr->pc.TDp;
        dTdH.at(g)= myPhase->pPtr->pc.dTD_dHD;
    }

    return 1;
}

// Effective velocity considering values on both sides of the edge.
// Interior edges
// Compute phase averaged velocity as well
int Driver::computeEffVel(const vector<vertex>& gaussp,
                          const vertexSet& edgep,
                          const indice& gcellin, const indice& gcellout,
                          const Tensor<weights>& allwgtsHD, double ** lHD,
                          const Tensor<weights>& allwgtsCD, double ** lCD,
                          vector<vertex>& effvel,
                          vector<vertex>& phasevel,
                          vector<vertex>& solidvel,
                          vector<double>& TDin,
                          vector<double>& TDout,
                          vector<double>& dTdHin,
                          vector<double>& dTdHout,
								  vector<double>& CDin,
								  vector<double>& CDout,
								  vector<double>& HDin,
								  vector<double>& HDout){

    // Extract porosity on both sides of the edges
    vector<double> phiin ; phiin.resize(gaussp.size());
    vector<double> phiout; phiout.resize(gaussp.size());
  
    // Extract concentration of solid and liquid on both sides of the edges
    vector<double> csout; csout.resize(gaussp.size());
    vector<double> clout; clout.resize(gaussp.size());
    vector<double> csin; csin.resize(gaussp.size());
    vector<double> clin; clin.resize(gaussp.size());

    // Extract velocity on these given gauss points
    vector<vertex> vel_relative = 
    ExtractVelocity(&sresult_->vel_darcy, &sresult_->g_darcy,
                    refArrayDarcyEssen_,mi,
                    gaussp, gcellin,*hdiv_,*basis_,{1});
    
    vector<vertex> vel_stokes = 
    ExtractVelocity(&sresult_->vel_stokes, &sresult_->g_stokes,
                    refArrayStokesEssen_,mi,
                    gaussp, gcellin,*br_,*basis_,{1});

    // Compute phase on both sides of the edge
    computephase(gaussp, edgep, gcellin, allwgtsHD, lHD, allwgtsCD, lCD,
                 csin, clin, phiin, TDin, dTdHin, CDin, HDin);

    computephase(gaussp, edgep, gcellout, allwgtsHD, lHD, allwgtsCD, lCD,
                 csout, clout, phiout, TDout, dTdHout, CDout, HDout);

    for (int g=0; g<gaussp.size(); g++){
        // Compute harmonic mean of physicalm values with respect to both sides 
        // of the edge 
        double cl_mean = 0.0;
        if (clin.at(g) == 0.0 && clout.at(g) == 0.0){
            cl_mean = 0.0;
        } else {
            cl_mean  = harmonic_mean(clin.at(g) ,clout.at(g));
        }

        double phi_mean = 0.0;
        if (phiin.at(g) == 0.0 && phiout.at(g) == 0.0){
            phi_mean = 0.0;
        } else {
            phi_mean  = harmonic_mean(phiin.at(g) ,phiout.at(g));
        }

        double cs_mean  = harmonic_mean(csin.at(g) ,csout.at(g));

        effvel.at(g) = cl_mean*phi_mean*(vel_relative.at(g) + vel_stokes.at(g)) + 
                       cs_mean*(1-phi_mean)*vel_stokes.at(g);
        effvel.at(g) /= cl_mean*phi_mean + cs_mean*(1-phi_mean);

        phasevel.at(g) = phi_mean*vel_relative.at(g) + vel_stokes.at(g);

        solidvel.at(g) = (1-phi_mean) * vel_stokes.at(g);

//        printf("vr %e, vs %e , phi %e , cl %e , cs %e , effvel %e , phasevel %e , solidvel %e\n", 
//              vel_relative.at(g)[1], vel_stokes.at(g)[1], phi_mean, cl_mean, 
//				  cs_mean, effvel.at(g)[1], phasevel.at(g)[1], solidvel.at(g)[1]);

        // Testing 
        //effvel.at(g)   = 1e-5;

        //phasevel.at(g) = 1e-5;
        //solidvel.at(g) = 1e-5;
    }

    return 1;
}

// One sided effective velocity
// Used on the boundary
int Driver::computeEffVel(const vector<vertex>& gaussp,
                          const vertexSet& edgep,
                          const indice& gcell,
                          const Tensor<weights>& allwgtsHD, double ** lHD,
                          const Tensor<weights>& allwgtsCD, double ** lCD,
                          vector<vertex>& effvel,
                          vector<vertex>& phasevel,
                          vector<vertex>& solidvel,
                          vector<double>& TD,
                          vector<double>& dTdH,
								  vector<double>& CD,
								  vector<double>& HD){

    // Extract porosity on both sides of the edges
    vector<double> phi ; phi.resize(gaussp.size());
  
    // Extract concentration of solid and liquid on both sides of the edges
    vector<double> cs; cs.resize(gaussp.size());
    vector<double> cl; cl.resize(gaussp.size());

    // Extract velocity on these given gauss points
    vector<vertex> vel_relative = 
    ExtractVelocity(&sresult_->vel_darcy, &sresult_->g_darcy,
                    refArrayDarcyEssen_,mi,
                    gaussp, gcell,*hdiv_,*basis_,{1});
    
    vector<vertex> vel_stokes = 
    ExtractVelocity(&sresult_->vel_stokes, &sresult_->g_stokes,
                    refArrayStokesEssen_,mi,
                    gaussp, gcell,*br_,*basis_,{1});

    // Compute phase on both sides of the edge
    computephase(gaussp, edgep, gcell, allwgtsHD, lHD, allwgtsCD, lCD,
                 cs, cl, phi, TD, dTdH, CD, HD);

    for (int g=0; g<gaussp.size(); g++){
        // Compute harmonic mean of physicalm values with respect to both sides 
        // of the edge 
        effvel.at(g) = cl.at(g)*phi.at(g)*(vel_relative.at(g) + vel_stokes.at(g)) + 
                       cs.at(g)*(1-phi.at(g))*vel_stokes.at(g);
        effvel.at(g) /= cl.at(g)*phi.at(g) + cs.at(g)*(1-phi.at(g));

        phasevel.at(g) = phi.at(g)*vel_relative.at(g) + vel_stokes.at(g);
        solidvel.at(g) = (1-phi.at(g)) * vel_stokes.at(g);

//        printf("vr %e, vs %e , phi %e , cl %e , cs %e \n", 
//              vel_relative.at(g)[1], vel_stokes.at(g)[1], phi.at(g), cl.at(g), 
//		   		  cs.at(g));

        // Testing 
        //effvel.at(g)   = 1e-5;

        //phasevel.at(g) = 1e-5;
        //solidvel.at(g) = 1e-5;
    } 

    return 1;
}

int Driver::updateEdgeFlux(Tensor<double>& vertedgeHD, Tensor<double>& horiedgeHD,
                           Tensor<double>& vertedgeCD, Tensor<double>& horiedgeCD,
                           const Tensor<weights>& allwgtsHD, double ** lHD,
                           const Tensor<weights>& allwgtsCD, double ** lCD){

    Tensor_zero(vertedgeCD);
    Tensor_zero(horiedgeCD);
    Tensor_zero(vertedgeHD);
    Tensor_zero(horiedgeHD);

    const valarray<double>& gwe = GaussWeightsEdge;
    const valarray<double>& gpe = GaussPointsEdge;

    std::vector<vertex> gaussp;
    gaussp.resize(gpe.size());

    // Effective velocity should be used for transport of concentration
    vector<vertex> effvel; effvel.resize(gaussp.size());
    vector<vertex> phasevel; phasevel.resize(gaussp.size());
    vector<vertex> solidvel; phasevel.resize(gaussp.size());
    vector<double> TDin; TDin.resize(gaussp.size());
    vector<double> TDout; TDout.resize(gaussp.size());
    vector<double> dTdHin; dTdHin.resize(gaussp.size());
    vector<double> dTdHout; dTdHout.resize(gaussp.size());

    vector<double> CDin; CDin.resize(gaussp.size());
    vector<double> CDout; CDout.resize(gaussp.size());
    vector<double> HDin; HDin.resize(gaussp.size());
    vector<double> HDout; HDout.resize(gaussp.size());
 
    for (int j=0; j<mi.MPIglobalCellSize[1]; j++){
    for (int i=0; i<mi.MPIglobalCellSize[0]; i++){

        double fluxHD = 0.0;
        double fluxCD = 0.0;
        double fluxL  = 0.0;

        indice gcell {i,j};
        indice cellout;

        clearall(effvel, phasevel, solidvel, TDin, TDout, dTdHin, dTdHout, 
                 CDin, CDout, HDin, HDout);

        // Extract corners with respect to given global indice
        vertexSet corners = extractCorners(mi, gcell); 

        // =============================================================
        // Get horizontal edge
        vertexSet hori {corners.at(0), corners.at(1)};
      
        // Extract velocity on this edge 
        for (int g=0; g<gpe.size(); g++){
            gaussp.at(g) = GaussMapPointsEdge({gpe[g]},hori);
        }   

        if(j==0){
            // Use one sided velocity
            computeEffVel(gaussp, hori, gcell, allwgtsHD, lHD, allwgtsCD, lCD, 
                          effvel, phasevel, solidvel, TDin, dTdHin, 
                          CDin, HDin);
            // bottom edge, fix values
            fluxCD = edgefluxintegral(hori, CDbottom, effvel);
            // =========================================================
            fluxHD = edgefluxintegral(hori, HDbottom, phasevel);
            fluxL  = myPhase->pPtr->LD * edgefluxintegral(hori, 1, solidvel);
        } else {
            // Use both sides velocity
            cellout = gcell + mi.faceNormal[0];
            computeEffVel(gaussp, hori, gcell, cellout, allwgtsHD, lHD, allwgtsCD, 
                          lCD, effvel, phasevel, solidvel, 
                          TDin, TDout, dTdHin, dTdHout, CDin, CDout, HDin, HDout);

//for (int g=0; g<gaussp.size(); g++){
//cout << effvel.at(g)[1] << "   ";
//}cout << endl;
            fluxCD = edgefluxintegral(mi, gcell, cellout, hori, allwgtsCD, 
                                      effvel, ml, advection, lCD);
//printf("%.16f, \n", fluxCD);
            // =========================================================
            fluxHD = edgefluxintegral(hori, HDin, HDout, TDin, TDout, 
                                      dTdHin, dTdHout, phasevel);
            fluxL  = myPhase->pPtr->LD * edgefluxintegral(hori, 1, solidvel);
        }
        horiedgeCD({i,j}) = fluxCD;
        horiedgeHD({i,j}) = fluxHD - fluxL;

        // =============================================================
        // Get vertical edge
        vertexSet vert {corners.at(3), corners.at(0)};

        gaussp.clear();gaussp.resize(gpe.size());
        // Extract velocity on this edge 
        for (int g=0; g<gpe.size(); g++){
            gaussp.at(g) = GaussMapPointsEdge({gpe[g]},vert);
        }   

        clearall(effvel, phasevel, solidvel, TDin, TDout, dTdHin, dTdHout, 
                 CDin, CDout, HDin, HDout);

        // boundary
        if (i==0){
            fluxCD = 0.0;
            fluxHD = 0.0;
            fluxL  = 0.0;
        } else {
            cellout = gcell + mi.faceNormal[3];
            computeEffVel(gaussp, vert, gcell, cellout, allwgtsHD, lHD, allwgtsCD, 
                          lCD, effvel, phasevel, solidvel, TDin, TDout, 
                          dTdHin, dTdHout, CDin, CDout, HDin, HDout);

            //fluxCD = edgefluxintegral(mi, gcell, cellout, vert, allwgtsCD, 
            //                          effvel, ml, advection, lCD);
            fluxHD = edgefluxintegral(vert, HDin, HDout, TDin, TDout, 
                                      dTdHin, dTdHout, phasevel);
            fluxL  = myPhase->pPtr->LD * edgefluxintegral(vert, 1, solidvel);
        }

        vertedgeCD({i,j}) = 0.0;//fluxCD;
        vertedgeHD({i,j}) = 0.0;//fluxHD - fluxL;

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

        clearall(effvel, phasevel, solidvel, TDin, TDout, dTdHin, dTdHout, 
                 CDin, CDout, HDin, HDout);

        for (int g=0; g<gpe.size(); g++){gaussp.at(g) = GaussMapPointsEdge({gpe[g]}, hori);}

        computeEffVel(gaussp, hori, gcell, allwgtsHD, lHD, allwgtsCD, lCD, 
                      effvel, phasevel, solidvel, TDin, dTdHin, CDin, HDin);

//        for (int g=0; g<gpe.size(); g++){printf("effvel %e, phasevel %e, solid vel %e ", effvel.at(g)[1], phasevel.at(g)[1], solidvel.at(g)[1]); cout << endl;}

        double fluxCD = edgefluxintegral(mi, gcell, hori, allwgtsCD, effvel, ml, advection, lCD);
        horiedgeCD({i, mi.MPIglobalCellSize[1]}) = fluxCD; 

        double fluxHD = edgefluxintegral(hori, HDin, HDin, TDin, TDin, 
                                       dTdHin, dTdHin, phasevel) - 
        myPhase->pPtr->LD * edgefluxintegral(hori, 1, solidvel);
  
        if (fluxHD > 0){
            fluxHD = 0.0;
        }

        horiedgeHD({i, mi.MPIglobalCellSize[1]}) = fluxHD; 
   }

    return 1;
}

int Driver::updateEdgeFlux(const Tensor<vertexSet>& phasevel_vert, 
                           const Tensor<vertexSet>& phasevel_hori, 
                           const Tensor<vertexSet>& effvel_vert, 
                           const Tensor<vertexSet>& effvel_hori, 
                           const Tensor<vertexSet>& solidvel_vert, 
                           const Tensor<vertexSet>& solidvel_hori,
                           Tensor<double>& vertedgeHD, Tensor<double>& horiedgeHD,
                           Tensor<double>& vertedgeCD, Tensor<double>& horiedgeCD,
                           const Tensor<weights>& allwgtsHD, double ** lHD,
                           const Tensor<weights>& allwgtsCD, double ** lCD){

    Tensor_zero(vertedgeCD);
    Tensor_zero(horiedgeCD);
    Tensor_zero(vertedgeHD);
    Tensor_zero(horiedgeHD);

    const valarray<double>& gwe = GaussWeightsEdge;
    const valarray<double>& gpe = GaussPointsEdge;

    std::vector<vertex> gaussp;
    gaussp.resize(gpe.size());

    vector<double> CDin; CDin.resize(gaussp.size());
    vector<double> CDout; CDout.resize(gaussp.size());
    vector<double> HDin; HDin.resize(gaussp.size());
    vector<double> HDout; HDout.resize(gaussp.size());
    vector<double> TDin; TDin.resize(gaussp.size());
    vector<double> TDout; TDout.resize(gaussp.size());
    vector<double> dTdHin; dTdHin.resize(gaussp.size());
    vector<double> dTdHout; dTdHout.resize(gaussp.size());

    vector<double> dummy; dummy.resize(gaussp.size());

    for (int j=0; j<mi.MPIglobalCellSize[1]; j++){
    for (int i=0; i<mi.MPIglobalCellSize[0]; i++){

        double fluxHD = 0.0;
        double fluxCD = 0.0;
        double fluxL  = 0.0;

        indice gcell {i,j};
        indice cellout;

        CDin.clear(); CDin.resize(gaussp.size());
        CDout.clear(); CDout.resize(gaussp.size());

        HDin.clear(); HDin.resize(gaussp.size());
        HDout.clear(); HDout.resize(gaussp.size());

        TDin.clear(); TDin.resize(gaussp.size());    
        TDout.clear(); TDout.resize(gaussp.size()); 

        dTdHin.clear(); dTdHin.resize(gaussp.size());    
        dTdHout.clear(); dTdHout.resize(gaussp.size()); 

        // Extract corners with respect to given global indice
        vertexSet corners = extractCorners(mi, gcell); 

        // =============================================================
        // Get horizontal edge
        vertexSet hori {corners.at(0), corners.at(1)};

        // Extract velocity on this edge 
        for (int g=0; g<gpe.size(); g++){
            gaussp.at(g) = GaussMapPointsEdge({gpe[g]},hori);
        }   

        computephase(gaussp, hori, gcell, allwgtsHD, lHD, allwgtsCD, lCD, 
                     dummy, dummy, dummy, TDin, dTdHin, CDin, HDin);

        if(j==0){
 
            // bottom edge, fix values
            fluxCD = edgefluxintegral(hori, CDbottom, effvel_hori({i,j}));
            // =========================================================
            fluxHD = edgefluxintegral(hori, HDbottom, phasevel_hori({i,j}));
            fluxL  = myPhase->pPtr->LD * 
                     edgefluxintegral(hori, 1, solidvel_hori({i,j}));
 
        } else {

            cellout = gcell + mi.faceNormal[0];

            computephase(gaussp, hori, cellout, allwgtsHD, lHD, allwgtsCD, lCD, 
                         dummy, dummy, dummy, TDout, dTdHout, CDout, HDout);
//for (int g=0; g<gaussp.size(); g++){
//cout << effvel_hori({i,j}).at(g)[1] << "   ";
//}cout << endl;
            fluxCD = edgefluxintegral(mi, gcell, cellout, hori, allwgtsCD, 
                                      effvel_hori({i,j}), ml, advection, lCD);
//printf("%.16f, \n", fluxCD);
            // =========================================================
            fluxHD = edgefluxintegral(hori, HDin, HDout, TDin, TDout, 
                                      dTdHin, dTdHout, phasevel_hori({i,j}));
            fluxL  = myPhase->pPtr->LD * 
                     edgefluxintegral(hori, 1, solidvel_hori({i,j}));

        }

        horiedgeCD({i,j}) = fluxCD;

        horiedgeHD({i,j}) = fluxHD - fluxL;

        // Cheating a little bit here

        vertedgeCD({i,j}) = 0.0;//fluxCD;
        vertedgeHD({i,j}) = 0.0;//fluxHD - fluxL;

    }}

    for (int i=0; i<mi.MPIglobalCellSize[0]; i++){

        // Regarded as outside cell
        indice gcell {i, mi.MPIglobalCellSize[1]-1};
        vertexSet corners = extractCorners(mi, gcell);
        vertexSet hori    = {corners.at(3), corners.at(2)};

        computephase(gaussp, hori, gcell, allwgtsHD, lHD, allwgtsCD, lCD, 
                     dummy, dummy, dummy, TDin, dTdHin, CDin, HDin);

        double fluxCD = edgefluxintegral(mi, gcell, hori, allwgtsCD, 
               effvel_hori({i,mi.MPIglobalCellSize[1]}), ml, advection, lCD);

        horiedgeCD({i, mi.MPIglobalCellSize[1]}) = fluxCD; 

        double fluxHD = edgefluxintegral(hori, HDin, HDin, TDin, TDin, 
               dTdHin, dTdHin, phasevel_hori({i,mi.MPIglobalCellSize[1]})) - 
        myPhase->pPtr->LD * edgefluxintegral(hori, 1, solidvel_hori({i,mi.MPIglobalCellSize[1]}));

        horiedgeHD({i, mi.MPIglobalCellSize[1]}) = fluxHD; 
    }

    return 1;
}

int Driver::computeFaceVel(const vector<vertex>& gaussp,
                           const indice& gcell, 
                           const Tensor<weights>& allwgtsHD, double ** lHD,
                           const Tensor<weights>& allwgtsCD, double ** lCD,
                           vector<vertex>& phasevel,
                           vector<double>& TD){

    double HD = 0.0, CD = 0.0;

    vector<vertex> vel_relative = 
    ExtractVelocity(&sresult_->vel_darcy, &sresult_->g_darcy,
                    refArrayDarcyEssen_,mi,
                    gaussp, gcell,*hdiv_,*basis_,{1});
    
    vector<vertex> vel_stokes = 
    ExtractVelocity(&sresult_->vel_stokes, &sresult_->g_stokes,
                    refArrayStokesEssen_,mi,
                    gaussp, gcell,*br_,*basis_,{1});

    for (int g=0; g<gaussp.size(); g++){

        HD = advection.eval(gaussp.at(g), ml, location(mi, gcell), 
             allwgtsHD({gcell[0], gcell[1]}), gcell, lHD);

        CD = advection.eval(gaussp.at(g), ml, location(mi, gcell), 
             allwgtsCD({gcell[0], gcell[1]}), gcell, lCD);

        double lithoP = myPhase->pPtr->GetStaticP(-1*gaussp.at(g)[1], 
                                                  myPhase->pPtr->l0); 
        myPhase->pPtr->evalPhase(HD, CD, lithoP);
 
        TD.at(g)  = myPhase->pPtr->pc.TDp;

        phasevel.at(g) = myPhase->pPtr->pc.phil * vel_relative.at(g) + 
                         vel_stokes.at(g);
    }

    return 1;
}

// Compute flux defined on face instead of edge
int Driver::updateCellFlux(Tensor<double>& fluxHD,
                           Tensor<double>& fluxCD,
                           const Tensor<weights>& allwgtsHD, double ** lHD,
                           const Tensor<weights>& allwgtsCD, double ** lCD){

    Tensor_zero(fluxHD);

    const valarray<double>& gwf = GaussWeightsFace;
    const vector<vertex>& gpf = GaussPointsFace;

    std::vector<vertex> gaussp;
    gaussp.resize(gwf.size());

    vector<vertex> phasevel; phasevel.resize(gaussp.size());
    vector<double> TD; TD.resize(gaussp.size());

    for (int j=0; j<mi.MPIglobalCellSize[1]; j++){
    for (int i=0; i<mi.MPIglobalCellSize[0]; i++){

        indice gcell {i,j};
        phasevel.clear(); phasevel.resize(gaussp.size());
        TD.clear(); TD.resize(gaussp.size());
        // Extract corners with respect to given global indice
        vertexSet corners = extractCorners(mi, gcell); 

        for (unsigned int g=0; g<gwf.size(); g++){
            gaussp.at(g) = GaussMapPointsFace(gpf[g],corners);
        }

        computeFaceVel(gaussp, gcell, allwgtsHD, lHD, allwgtsCD, lCD, phasevel, TD);

        double work = 0.0;

        for (unsigned int g=0; g<gwf.size(); g++){

            double jac = abs(GaussJacobian(gpf[g],corners));
            double gw = gwf[g];
    
            work += -10*phasevel.at(g)[1]*TD.at(g)* jac * gw;

        }

        double area = mi.cellArea.at(FlatIndic(mi,gcell));

        fluxHD({i,j}) = work*myPhase->pPtr->alpha0 * myPhase->pPtr->l0 / 
                        myPhase->pPtr->cp / area;

    }}

    return 1;
}

int Driver::updateVel_Pause(Tensor<vertexSet>& phasevel_vert, 
                            Tensor<vertexSet>& phasevel_hori, 
                            Tensor<vertexSet>& effvel_vert, 
                            Tensor<vertexSet>& effvel_hori, 
                            Tensor<vertexSet>& solidvel_vert, 
                            Tensor<vertexSet>& solidvel_hori, 
                            const Tensor<weights>& allwgtsHD, double ** lHD,
                            const Tensor<weights>& allwgtsCD, double ** lCD){

    const valarray<double>& gwe = GaussWeightsEdge;
    const valarray<double>& gpe = GaussPointsEdge;

    std::vector<vertex> gaussp;
    gaussp.resize(gpe.size());

    // Effective velocity should be used for transport of concentration
    vector<vertex> effvel; effvel.resize(gaussp.size());
    vector<vertex> phasevel; phasevel.resize(gaussp.size());
    vector<vertex> solidvel; phasevel.resize(gaussp.size());
    vector<double> TDin; TDin.resize(gaussp.size());
    vector<double> TDout; TDout.resize(gaussp.size());
    vector<double> dTdHin; dTdHin.resize(gaussp.size());
    vector<double> dTdHout; dTdHout.resize(gaussp.size());

    vector<double> CDin; CDin.resize(gaussp.size());
    vector<double> CDout; CDout.resize(gaussp.size());
    vector<double> HDin; HDin.resize(gaussp.size());
    vector<double> HDout; HDout.resize(gaussp.size());

    for (int j=0; j<mi.MPIglobalCellSize[1]; j++){
    for (int i=0; i<mi.MPIglobalCellSize[0]; i++){

        indice gcell {i,j};
        indice cellout;

        clearall(effvel, phasevel, solidvel, TDin, TDout, dTdHin, dTdHout, 
                 CDin, CDout, HDin, HDout);
 
        // Extract corners with respect to given global indice
        vertexSet corners = extractCorners(mi, gcell); 

        // =============================================================
        // Get horizontal edge
        vertexSet hori {corners.at(0), corners.at(1)};
      
        // Extract velocity on this edge 
        for (int g=0; g<gpe.size(); g++){
            gaussp.at(g) = GaussMapPointsEdge({gpe[g]},hori);
        }   

        if(j==0){
            // Use one sided velocity
            computeEffVel(gaussp, hori, gcell, allwgtsHD, lHD, allwgtsCD, lCD, 
                          effvel, phasevel, solidvel, TDin, dTdHin, 
                          CDin, HDin);
        } else {

            cellout = gcell + mi.faceNormal[0];
            computeEffVel(gaussp, hori, gcell, cellout, allwgtsHD, lHD, allwgtsCD, 
                          lCD, effvel, phasevel, solidvel, 
                          TDin, TDout, dTdHin, dTdHout, CDin, CDout, HDin, HDout);
        }

            phasevel_hori({i,j}) = phasevel;
            effvel_hori({i,j}) = effvel;
            solidvel_hori({i,j}) = solidvel;

        // =============================================================
        // Get vertical edge
        vertexSet vert {corners.at(3), corners.at(0)};

        gaussp.clear();gaussp.resize(gpe.size());
        // Extract velocity on this edge 
        for (int g=0; g<gpe.size(); g++){
            gaussp.at(g) = GaussMapPointsEdge({gpe[g]},vert);
        }   

        clearall(effvel, phasevel, solidvel, TDin, TDout, dTdHin, dTdHout, 
                 CDin, CDout, HDin, HDout);

        if (i==0){
            // Use one sided velocity
            computeEffVel(gaussp, vert, gcell, allwgtsHD, lHD, allwgtsCD, lCD, 
                          effvel, phasevel, solidvel, TDin, dTdHin, 
                          CDin, HDin);

        } else {
            cellout = gcell + mi.faceNormal[3];
            computeEffVel(gaussp, vert, gcell, cellout, allwgtsHD, lHD, allwgtsCD, 
                          lCD, effvel, phasevel, solidvel, TDin, TDout, 
                          dTdHin, dTdHout, CDin, CDout, HDin, HDout);
        }

        phasevel_vert({i,j}) = phasevel;
        effvel_vert({i,j}) = effvel;
        solidvel_vert({i,j}) = solidvel;

    }}

    for (int j=0; j<mi.MPIglobalCellSize[1]; j++){

        indice gcell {mi.MPIglobalCellSize[0]-1, j};
        vertexSet corners = extractCorners(mi, gcell);
        vertexSet vert    = {corners.at(2), corners.at(1)};

        clearall(effvel, phasevel, solidvel, TDin, TDout, dTdHin, dTdHout, 
                 CDin, CDout, HDin, HDout);

        for (int g=0; g<gpe.size(); g++){gaussp.at(g) = GaussMapPointsEdge({gpe[g]}, vert);}

        computeEffVel(gaussp, vert, gcell, allwgtsHD, lHD, allwgtsCD, lCD, 
                      effvel, phasevel, solidvel, TDin, dTdHin, CDin, HDin);

        phasevel_vert({mi.MPIglobalCellSize[0],j}) = phasevel;
        effvel_vert({mi.MPIglobalCellSize[0],j}) = effvel;
        solidvel_vert({mi.MPIglobalCellSize[0],j}) = solidvel;
    }

    // right and top edges
    for (int i=0; i<mi.MPIglobalCellSize[0]; i++){
        // Regarded as outside cell
        indice gcell {i, mi.MPIglobalCellSize[1]-1};
        vertexSet corners = extractCorners(mi, gcell);
        vertexSet hori    = {corners.at(3), corners.at(2)};

        clearall(effvel, phasevel, solidvel, TDin, TDout, dTdHin, dTdHout, 
                 CDin, CDout, HDin, HDout);

        for (int g=0; g<gpe.size(); g++){gaussp.at(g) = GaussMapPointsEdge({gpe[g]}, hori);}

        computeEffVel(gaussp, hori, gcell, allwgtsHD, lHD, allwgtsCD, lCD, 
                      effvel, phasevel, solidvel, TDin, dTdHin, CDin, HDin);

        phasevel_hori({i,mi.MPIglobalCellSize[1]}) = phasevel;
        effvel_hori({i,mi.MPIglobalCellSize[1]}) = effvel;
        solidvel_hori({i,mi.MPIglobalCellSize[1]}) = solidvel;
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

    Tensor<double> facefluxHD = Tensor<double>(2);
    facefluxHD.setSize({mi.MPIlocalCellSize[0], mi.MPIlocalCellSize[1]});
    Tensor<double> facefluxCD = Tensor<double>(2);
    facefluxCD.setSize({mi.MPIlocalCellSize[0], mi.MPIlocalCellSize[1]});

    updateCellFlux(facefluxHD, facefluxCD, allwgtsHD, lHD, allwgtsCD, lCD);

    for (int j=0; j<mi.MPIglobalCellSize[1]; j++){
    for (int i=0; i<mi.MPIglobalCellSize[0]; i++){

        lfHD[j][i] = getcellflux(mi, {i,j}, vertedgefluxHD, horiedgefluxHD) ;
					 //- facefluxHD({i,j});
        lfCD[j][i] = getcellflux(mi, {i,j}, vertedgefluxCD, horiedgefluxCD);

//        cout << lfHD[j][i] << "   ";
//        cout << lfCD[j][i] << "   ";

    }}
    return 1;
}

int Driver::getflux(const Tensor<vertexSet>& phasevel_vert, 
                    const Tensor<vertexSet>& phasevel_hori, 
                    const Tensor<vertexSet>& effvel_vert, 
                    const Tensor<vertexSet>& effvel_hori, 
                    const Tensor<vertexSet>& solidvel_vert, 
                    const Tensor<vertexSet>& solidvel_hori,
                    const Tensor<weights>& allwgtsHD, double ** lHD,
                    const Tensor<weights>& allwgtsCD, double ** lCD,
                    double **lfHD, double **lfCD){

    Tensor<double> horiedgefluxHD = Tensor<double>(2);
    horiedgefluxHD.setSize({mi.MPIlocalCellSize[0], mi.MPIlocalCellSize[1]+1});

    Tensor<double> vertedgefluxHD = Tensor<double>(2);
    vertedgefluxHD.setSize({mi.MPIlocalCellSize[0]+1, mi.MPIlocalCellSize[1]});

    Tensor<double> horiedgefluxCD = Tensor<double>(2);
    horiedgefluxCD.setSize({mi.MPIlocalCellSize[0], mi.MPIlocalCellSize[1]+1});

    Tensor<double> vertedgefluxCD = Tensor<double>(2);
    vertedgefluxCD.setSize({mi.MPIlocalCellSize[0]+1, mi.MPIlocalCellSize[1]});

    updateEdgeFlux(phasevel_vert, phasevel_hori, 
                   effvel_vert,   effvel_hori,
                   solidvel_vert, solidvel_hori,
                   vertedgefluxHD, horiedgefluxHD,
                   vertedgefluxCD, horiedgefluxCD,
                   allwgtsHD, lHD,
                   allwgtsCD, lCD);

    //Tensor<double> facefluxHD = Tensor<double>(2);
    //facefluxHD.setSize({mi.MPIlocalCellSize[0], mi.MPIlocalCellSize[1]});
    //Tensor<double> facefluxCD = Tensor<double>(2);
    //facefluxCD.setSize({mi.MPIlocalCellSize[0], mi.MPIlocalCellSize[1]});

    //updateCellFlux(facefluxHD, facefluxCD, allwgtsHD, lHD, allwgtsCD, lCD);

    for (int j=0; j<mi.MPIglobalCellSize[1]; j++){
    for (int i=0; i<mi.MPIglobalCellSize[0]; i++){

        lfHD[j][i] = getcellflux(mi, {i,j}, vertedgefluxHD, horiedgefluxHD);
					 //- facefluxHD({i,j});
        lfCD[j][i] = getcellflux(mi, {i,j}, vertedgefluxCD, horiedgefluxCD);

        //cout << lfHD[j][i] << "   " ;
        //cout << lfCD[j][i] << "   " ;
    } }
 

    return 1;
}
