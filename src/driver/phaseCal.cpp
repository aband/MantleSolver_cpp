#include "phaseCal.h"

inline double evalCbar2(const double& CD,
                        const double& HD,
                        const double& cbar,
                        EUTECTIC::evalPhase* pPtr){

    return  pPtr->getD1s()*pPtr->getPhi2() + 
            pPtr->getD2f()/pPtr->getD1f()  * 
           (pPtr->getD1f()*pPtr->getPhif() + 
            pPtr->getD1s()*pPtr->getPhi1() - cbar); 
}

inline double evalDcbar2(const double& CD,
                         const double& HD,
                         const double& cbar,
                         EUTECTIC::evalPhase* pPtr){

    return pPtr->getD1s()*pPtr->getDPhi2() + 
           pPtr->getD2f()/pPtr->getD1f()  * 
          (pPtr->getD1f()*pPtr->getDPhif() + 
           pPtr->getD1s()*pPtr->getDPhi1() - cbar); 
}

inline double evalFunc(const double& CD,
                       const double& HD,
                       const double& cbar,
                       EUTECTIC::evalPhase* pPtr){

    return CD - 1 + cbar / (cbar + evalCbar2(CD,HD,cbar,pPtr));
}

inline double evalDFunc(const double& CD,
                        const double& HD,
                        const double& cbar,
                        EUTECTIC::evalPhase* pPtr){

    return 1 - cbar*evalDcbar2(CD,HD,cbar,pPtr)/ pow(cbar+evalCbar2(CD,HD,cbar,pPtr),2);

}

double ComputePhase(const double& initCD, 
                    const double& HD,
                    const double& cbar,
                    EUTECTIC::evalPhase* pPtr){

    // Compute current phase consistent to phase package
    // Evaluate with evalPhase pointer 

    // Evaluate current phase behavior first
    pPtr->EvalPhase(HD,initCD);

    double nextCD = initCD;
    double oldCD = initCD;
    double MAX = 10;
    double tol = 1e-5;

    for (int count = 0; count < MAX; count ++){
        nextCD = nextCD - evalFunc(nextCD,HD,cbar,pPtr)/evalDFunc(nextCD,HD,cbar,pPtr);
		  std::cout << nextCD - oldCD << std::endl;
        //if (abs(nextCD - oldCD) < tol){
		//			 std::cout << count << std::endl;
       // }
        oldCD = nextCD;
    }

    return nextCD;
}
