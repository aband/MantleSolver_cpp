#include "driver.h"

int clearall(vector<vertex>& effvel,
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

// Compute phase values on a given edge with respect to 
// a given cell
int Driver::computephase(const std::vector<vertex>& gaussp,
                         const vertexSet& edgep,
                         const indice& gcell,
                         double ** lHD, double ** lCD,
                         vector<double>& cs,
                         vector<double>& cl,
                         vector<double>& phi,
                         vector<double>& TD,
                         vector<double>& dTdH,
	                      vector<double>& cd,
                         vector<double>& hd){

    double HD = 0.0, CD = 0.0;

    int s = FlatIndic(mi, gcell);

    for (int g=0; g<gaussp.size(); g++){

        HD = my_recon_HD.at(s).eval(lHD, gaussp.at(g), stenlg, stensm);
        CD = my_recon_CD.at(s).eval(lCD, gaussp.at(g), stenlg, stensm);

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
                          double ** lHD, double ** lCD,
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
    computephase(gaussp, edgep, gcellin, lHD, lCD,
                 csin, clin, phiin, TDin, dTdHin, CDin, HDin);

    computephase(gaussp, edgep, gcellout, lHD, lCD,
                 csout, clout, phiout, TDout, dTdHout, CDout, HDout);

    for (int g=0; g<gaussp.size(); g++){
        // Compute harmonic mean of physicalm values with respect to both sides 
        // of the edge 
        double cl_mean = 0.0;
        if (clin.at(g) == 0.0 && clout.at(g) == 0.0){
            cl_mean = 0.0;
        } else {
            //cl_mean  = harmonic_mean(clin.at(g) ,clout.at(g));
            cl_mean = (clin.at(g) + clout.at(g))/2.0;
        }

        double phi_mean = 0.0;
        if (phiin.at(g) == 0.0 && phiout.at(g) == 0.0){
            phi_mean = 0.0;
        } else {
            phi_mean  = harmonic_mean(phiin.at(g) ,phiout.at(g));
        }

        double cs_mean  = harmonic_mean(csin.at(g) ,csout.at(g));
        cs_mean  = (csin.at(g) +csout.at(g))/2.0;

        effvel.at(g) = cl_mean*phi_mean*(vel_relative.at(g) + vel_stokes.at(g)) + 
                       cs_mean*(1-phi_mean)*vel_stokes.at(g);
        effvel.at(g) /= cl_mean*phi_mean + cs_mean*(1-phi_mean);
//cout << effvel.at(g)[0] << "  " << effvel.at(g)[1] << endl;
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
                          double ** lHD,
                          double ** lCD,
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
    computephase(gaussp, edgep, gcell, lHD, lCD,
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

// Effective velocity considering values on both sides of the edge.
// Interior edges
// Compute phase averaged velocity as well
int Driver::computeEffVel_Nonlinear(const vector<vertex>& gaussp,
                                    const vertexSet& edgep,
                                    const indice& gcellin, const indice& gcellout,
                                    double ** lHD, double ** lCD,
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
			          					   vector<double>& HDout,
											   vector<double>& nonlinuin,
											   vector<double>& nonlinuout,
											   vector<double>& nonlinfin,
											   vector<double>& nonlinfout,
											   vector<double>& nonlindfduin,
												vector<double>& nonlindfduout,
												vector<vertex>& nonlinvel){

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
    computephase(gaussp, edgep, gcellin, lHD, lCD,
                 csin, clin, phiin, TDin, dTdHin, CDin, HDin);

    computephase(gaussp, edgep, gcellout, lHD, lCD,
                 csout, clout, phiout, TDout, dTdHout, CDout, HDout);

    for (int g=0; g<gaussp.size(); g++){
        // Compute harmonic mean of physicalm values with respect to both sides 
        // of the edge 
        double cl_mean = 0.0;
        if (clin.at(g) == 0.0 && clout.at(g) == 0.0){
            cl_mean = 0.0;
        } else {
            //cl_mean  = harmonic_mean(clin.at(g) ,clout.at(g));
            cl_mean = (clin.at(g) + clout.at(g))/2.0;
        }

        double phi_mean = 0.0;
        if (phiin.at(g) == 0.0 && phiout.at(g) == 0.0){
            phi_mean = 0.0;
        } else {
            phi_mean  = harmonic_mean(phiin.at(g) ,phiout.at(g));
        }

        double cs_mean  = harmonic_mean(csin.at(g) ,csout.at(g));
        cs_mean  = (csin.at(g) +csout.at(g))/2.0;

        effvel.at(g) = cl_mean*phi_mean*(vel_relative.at(g) + vel_stokes.at(g)) + 
                       cs_mean*(1-phi_mean)*vel_stokes.at(g);
        effvel.at(g) /= cl_mean*phi_mean + cs_mean*(1-phi_mean);
//cout << effvel.at(g)[0] << "  " << effvel.at(g)[1] << endl;
        phasevel.at(g) = phi_mean*vel_relative.at(g) + vel_stokes.at(g);

        solidvel.at(g) = (1-phi_mean) * vel_stokes.at(g);

//        printf("vr %e, vs %e , phi %e , cl %e , cs %e , effvel %e , phasevel %e , solidvel %e\n", 
//              vel_relative.at(g)[1], vel_stokes.at(g)[1], phi_mean, cl_mean, 
//				  cs_mean, effvel.at(g)[1], phasevel.at(g)[1], solidvel.at(g)[1]);

        // Testing 
        //effvel.at(g)   = 1e-5;

        //phasevel.at(g) = 1e-5;
        //solidvel.at(g) = 1e-5;

        // ================= Get nonlinear flux at the same time
        vertex tmp = clin.at(g)*phi_mean*(vel_relative.at(g) + vel_stokes.at(g))  + 
                     csin.at(g)*(1-phi_mean)*vel_stokes.at(g);
        tmp /= clin.at(g)*phi_mean + csin.at(g)*(1-phi_mean);

        nonlinuin.at(g) = CDin.at(g);
        nonlinfin.at(g) = tmp[1]*CDin.at(g);
   
        nonlindfduin.at(g) = tmp[1] + 1e-5;
  
        tmp = clout.at(g)*phi_mean*(vel_relative.at(g) + vel_stokes.at(g))  + 
              csout.at(g)*(1-phi_mean)*vel_stokes.at(g);
        tmp /= clout.at(g)*phi_mean + csout.at(g)*(1-phi_mean);

        nonlinuout.at(g) = CDout.at(g);
        nonlinfout.at(g) = tmp[1]*CDout.at(g);

        nonlindfduout.at(g) = tmp[1] + 1e-5;        

        nonlinvel.at(g) = {0,1};
    }

    return 1;
}

int Driver::getfluxall(Vec * fHDin, Vec * fCDin, bool updateVel, double t){

    int M = mi.MPIglobalCellSize[0];
    int N = mi.MPIglobalCellSize[1];
  
    vector<double> sigma_lg_HD;
    sigma_lg_HD.resize(stenlg.size());

    vector<double> sigma_sm_HD;
    sigma_sm_HD.resize(stensm.size());

    vector<double> sigma_lg_CD;
    sigma_lg_CD.resize(stenlg.size());

    vector<double> sigma_sm_CD;
    sigma_sm_CD.resize(stensm.size());

    Vec localHD, localCD; 
    PetscCall(DMGetLocalVector(dmu, &localHD));
    PetscCall(DMGetLocalVector(dmu, &localCD));

    PetscCall(DMGlobalToLocalBegin(dmu, globalHD, INSERT_VALUES, localHD));
    PetscCall(DMGlobalToLocalEnd(dmu, globalHD, INSERT_VALUES, localHD));

    PetscCall(DMGlobalToLocalBegin(dmu, globalCD, INSERT_VALUES, localCD));
    PetscCall(DMGlobalToLocalEnd(dmu, globalCD, INSERT_VALUES, localCD));

    Vec fluxHD = *fHDin;
    Vec fluxCD = *fCDin;

    double ** lHD;
    DMDAVecGetArray(dmu, localHD, &lHD);

    double ** fHD;
    DMDAVecGetArray(dmu, fluxHD, &fHD);

    double ** lCD;
    DMDAVecGetArray(dmu, localCD, &lCD);

    double ** fCD;
    DMDAVecGetArray(dmu, fluxCD, &fCD);

    for (int s=0; s<stenlg.size(); s++){
        sigma_lg_HD.at(s) = stenlg.at(s).sigma(lHD);
        sigma_lg_CD.at(s) = stenlg.at(s).sigma(lCD);
    }

    for (int s=0; s<stensm.size(); s++){
        sigma_sm_HD.at(s) = stensm.at(s).sigma(lHD);
        sigma_sm_CD.at(s) = stensm.at(s).sigma(lCD);
    }

    // Setup nonlinear weights
    for (int s=0; s<my_recon_HD.size(); s++){
        my_recon_CD.at(s).extractsigma(sigma_lg_CD, sigma_sm_CD);
        my_recon_CD.at(s).setWgts(L_*H_/(double)M/(double)N);

        my_recon_HD.at(s).extractsigma(sigma_lg_HD, sigma_sm_HD);
        my_recon_HD.at(s).setWgts(L_*H_/(double)M/(double)N);
    }

    if (updateVel){
        SolveFlow(maxIter, tolUzawa, lHD, lCD);
        CreateScatterVec();
    }

    vector<double> edgefluxHD;
    vector<double> edgefluxCD;

    computeEdgeFlux(edgefluxHD, edgefluxCD, t, lHD, lCD);

    for (int j=0; j<N; j++){
    for (int i=0; i<M; i++){

       double area = mi.cellArea.at(FlatIndic(mi, {i,j}));

       int bottom = j*M+i;
       int top    = (j+1)*M+i;
       int left   = M*(N+1) + j*(M+1) + i;
       int right  = M*(N+1) + j*(M+1) + i+1;

       fHD[j][i] = (edgefluxHD.at(bottom) - edgefluxHD.at(top) + 
                    edgefluxHD.at(left) - edgefluxHD.at(right))/area;

       fCD[j][i] = (edgefluxCD.at(bottom) - edgefluxCD.at(top) + 
                    edgefluxCD.at(left) - edgefluxCD.at(right))/area;
    }}

    DMDAVecRestoreArray(dmu, fluxHD, &fHD);
    DMDAVecRestoreArray(dmu, fluxCD, &fCD);

    DMDAVecRestoreArray(dmu, localHD, &lHD);
    DMRestoreLocalVector(dmu, &localHD);

    DMDAVecRestoreArray(dmu, localCD, &lCD);
    DMRestoreLocalVector(dmu, &localCD);

    return 1;
}
