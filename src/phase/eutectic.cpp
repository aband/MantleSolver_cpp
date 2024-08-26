phaseEval::phase(double nP){

    // Set up Clapeyron constant
    gamma_ = ;

    // Set up melting points under standard atmospheric pressure
	 // with dimension
	 te0_ = 1560; //(K) 
	 t10_ = 2053; //(K) 

	 // dimensionless
	 // regardless of value of pressure
    Te_ = 0;
	 T1_ = 1;
	
	 // dimensionless latent heat 
    L_ = 0.5;

}

void evalPhase(){


}

int phaseEval::phaseSplit(const double& H, 
                          const double& C){


    return 0;
}
