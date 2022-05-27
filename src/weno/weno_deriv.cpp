void WenoReconst::CreateOmegaHatDerivative(const WenoMesh*& wm){

    omega_deriv_lwp = new double[Lorder[0]*Lorder[1]];

    for (int s=0; s < omega_deriv_swp.size(); s++){
        omega_deriv_swp[s] = new double[Lorder[0]*Lorder[1]];
    }

    for (int e=0; e<StencilLarge.size(); e++){

    }

}

double * WenoReconst::CreateOmegaTildeDerivative(const WenoMesh*& wm, point target_point){

    double * domegat = new double[Lorder[0]*Lorder[1]];

    return domeget; 
}

void WenoReconst::CreateSmoothnessIndDerivative(const WenoMesh*& wm, point target_point){



}

double * WenoReconst::CreateWenoDerivative(const WenoMesh*& wm, point target_point){

    double * wd = new double[Lorder[0]*Lorder[1]];

    // Define offset parameters
    int offseti = StencilLarge[0][0];
    int offsetj = StencilLarge[0][1];

    int N = Lorder[0];

    int n = Lorder[0]*Lorder[1];
    double phi = 0.0;

    for (int p=0; p<n; p++){
        int locali = StencilLarge[p][0] - offseti;
        int localj = StencilLarge[p][1] - offsetj;

        int target_i = target[0] - wm->ghost + StencilLarge[p][0];

    }

    return wd;
}
