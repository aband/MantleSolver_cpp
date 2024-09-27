#include "driver.h"

int Driver::PrepareFlow(){

    basis_ = new basis();
    hdiv_  = new Hdivmixed();
    br_    = new BRMixed();

    br_->ComputeTotalDOF(mi);
    hdiv_->ComputeTotalDOF(mi);

    MarkBndryDOFStokes(bndryStokesEssen_, bndryStokesNatur_, mi, *basis_, *br_, myPhase->pp);
    MarkBndryDOFDarcy(bndryDarcyEssen_, bndryDarcyNatur_, mi, *basis_, *hdiv_, myPhase->pp);

    reducedDarcy_ = (ReducedSys *)malloc(sizeof(ReducedSys));
    reducedStokes_ = (ReducedSys *)malloc(sizeof(ReducedSys));

    refArrayStokes_ = new int[br_->getDOF()];
    refArrayDarcy_  = new int[hdiv_->getDOF()];

    return 1;
}

int Driver::SolveFlow(){

    ParallelMatrixAssemble(mi, *basis_, myPhase, bndryStokesEssen_, reducedStokes_, 
                                                 bndryDarcyEssen_,  reducedDarcy_, 
                           &K_, *br_, *hdiv_ , mluseAdv_,

                           refArrayStokes_, refArrayDarcy_, bndryDOFStokes_, bndryDOFDarcy_);

    return 1;
}
