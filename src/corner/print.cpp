#include "couple.h"

int couple::printGaussPoints(){

    FILE * vertggridx = fopen("vertgaussgridx.dat", "w");
    FILE * vertggridy = fopen("vertgaussgridy.dat", "w");

    FILE * horiggridx = fopen("horigaussgridx.dat", "w");
    FILE * horiggridy = fopen("horigaussgridy.dat", "w");

    // Vertical points first
    for (int j=0; j<N_  ; j++){
    for (int i=0; i<M_+1; i++){
      
        int dof = (j*(M_+1) + i)*3;
 
        for(int g=0; g<3; g++){
            fprintf(vertggridx, "%e ", edgegauss.at(dof+g)[0]);
            fprintf(vertggridy, "%e ", edgegauss.at(dof+g)[1]);
				//cout <<dof+g << "  " << edgegauss.at(dof + g)[0] << "  " << edgegauss.at(dof + g )[1] << endl;
        }
    }}
//    }fprintf(vertggridx, "\n ");
//     fprintf(vertggridy, "\n ");}

    int tolvert = N_*(M_+1)*3;

    // Horizontal points second
    for (int j=0; j<N_+1; j++){
    for (int i=0; i<M_;   i++){

        int dof = tolvert + (j*M_ + i)*3;
        for (int g=0; g<3; g++){
            fprintf(horiggridx, "%e ", edgegauss.at(dof+g)[0]);
            fprintf(horiggridy, "%e ", edgegauss.at(dof+g)[1]);
        }
    }}

//    }fprintf(horiggridx, "\n ");
//     fprintf(horiggridy, "\n ");}

    fclose(vertggridx);
    fclose(vertggridy);
    fclose(horiggridx);
    fclose(horiggridy);

    return 1;
}

int couple::printCellGrids(){

    FILE * cellgridx = fopen("cellgridx.dat", "w");
    FILE * cellgridy = fopen("cellgridy.dat", "w");

    for (int j=0; j<N_; j++){
    for (int i=0; i<M_; i++){

        fprintf(cellgridx, "%e ", cellcenter.at(j*M_+i)[0]);
        fprintf(cellgridy, "%e ", cellcenter.at(j*M_+i)[1]);

    }}

    fclose(cellgridx);
    fclose(cellgridy);
    return 1;
}

char * GetFilename(const char * fieldname, int mark){

    char * filename = (char *)malloc(strlen(fieldname)+10+4);

    char n_char[10];
    std::sprintf(n_char,"%d",mark);
    strcpy(filename, fieldname);
    strcat(filename, n_char);
    strcat(filename, ".dat");

    return filename;
}

char * GetFilenameAdd(const char * fieldname, const char * add, int mark){

    char * filename = (char *)malloc(strlen(fieldname)+15+4);

    char n_char[15];
    std::sprintf(n_char,"%d",mark);
    strcpy(filename, fieldname);
    strcat(filename, add);
    strcat(filename, n_char);
    strcat(filename, ".dat");

    return filename;
}

int couple::printCellScalar(Vec * sol, const char * fieldname, int mark){

    Vec temp = *sol;

    FILE * file = fopen(GetFilename(fieldname, mark), "w");

    for(int j=0; j<N_; j++){
    for(int i=0; i<M_; i++){

        double val;
        indice index {i,j};

        int nelem = FlatIndic(mi, index);

        PetscCall(VecGetValues(temp,1, &nelem, &val));

        fprintf(file, "%.14lf ", val);

    }fprintf(file, "\n");}

    fclose(file);

    return 1;
}

int couple::printedgeval(int mark, const vector<double>& val,
                                   const char * fieldname){

    FILE * file = fopen(GetFilename(fieldname, mark), "w");

    // Print scalar values on gauss quadrature points
	 // Vertical points first
    for (int j=0; j<N_  ; j++){
    for (int i=0; i<M_+1; i++){

        int dof = (j*(M_+1) + i)*3;

        for (int g=0; g<3; g++){
            fprintf(file,"%12f ", val.at(dof+g));
        }

    }}

    int tolvert = N_*(M_+1)*3;

    // Horizontal points next
    for (int j=0; j<N_+1; j++){
    for (int i=0; i<M_;   i++){

        int dof = tolvert + (j*M_ + i)*3;
        for (int g=0; g<3; g++){
            fprintf(file, "%12f ", val.at(dof+g));
        }
    }}

    fclose(file);
    return 1;
}

int couple::printcellval(int mark, const vector<double>& val,
                                   const char * fieldname){

    // Print values assigned to cell center
    FILE * file = fopen(GetFilename(fieldname, mark), "w");

    for (int j=0; j<N_; j++){
    for (int i=0; i<M_; i++){

        fprintf(file, "%12f ", val.at(j*M_+i));

    }}

    fclose(file);
    return 1;
}

int couple::printedgeval(int mark, const vector<vertex>& val,
                                   const char * fieldname){

    FILE * filex = fopen(GetFilenameAdd(fieldname, "x", mark), "w");
    FILE * filey = fopen(GetFilenameAdd(fieldname, "y", mark), "w");

    // Print scalar values on gauss quadrature points
    // Vertical points first
    for (int j=0; j<N_  ; j++){
    for (int i=0; i<M_+1; i++){

        int dof = (j*(M_+1) + i)*3;

        for (int g=0; g<3; g++){
            fprintf(filex,"%12f ", val.at(dof+g)[0]);
            fprintf(filey,"%12f ", val.at(dof+g)[1]);
        }

    }}

    int tolvert = N_*(M_+1)*3;

    // Horizontal points next
    for (int j=0; j<N_+1; j++){
    for (int i=0; i<M_;   i++){

        int dof = tolvert + (j*M_ + i)*3;
        for (int g=0; g<3; g++){
            fprintf(filex, "%12f ", val.at(dof+g)[0]);
            fprintf(filey, "%12f ", val.at(dof+g)[1]);
        }
    }}

    fclose(filex);
    fclose(filey);

    return 1;
}

int couple::printedgeporosity(int mark){
cout << "here" << endl;
    printedgeval(mark, edgeporo, "porosity");

    printcellval(mark, average_poro, "aveporo");

    return 1;
}

int couple::printphase(int mark){

    FILE * file = fopen(GetFilename("phase", mark), "w");
    FILE * filetdp = fopen(GetFilename("currentp", mark), "w");
    FILE * filetmp = fopen(GetFilename("meltp", mark), "w");

    for (int j=0; j<N_; j++){
    for (int i=0; i<M_; i++){

        fprintf(file, "%d ", cellphase.at(j*M_+i));
        fprintf(filetdp, "%lf ", currentTemp.at(j*M_+i));
        fprintf(filetmp, "%lf ", meltp.at(j*M_+i));
	 }}

    fclose(file);
    fclose(filetdp);
    fclose(filetmp);

    return 1;
}

int couple::printedgevel(int mark, const vector<vertex>& stokesvel,
                                   const vector<vertex>& darcyvel){

    printedgeval(mark, stokesvel, "edgevel_stokes"); 

    printedgeval(mark, darcyvel , "edgevel_darcy");

    // Print porosity corrected darcy velocity
    vector<vertex> target;
    target.resize(darcyvel.size());

    for (int s=0; s<(int)darcyvel.size(); s++){
        target.at(s) = darcyvel.at(s)*edgeporo.at(s);
	 }

    printedgeval(mark, target, "edgevel_darcy_porosity");

    for (int s=0; s<(int)darcyvel.size(); s++){
        target.at(s) = target.at(s) + stokesvel.at(s);
	 }

    printedgeval(mark, target, "phaseave_vel");

    return 1;
}

int couple::printpressure(int mark, Vec * sp, Vec * dp, bool printpressure, 
					 const char * spname, const char * dpname){

    // Break pressure vector

    Vec vecq       = *sp;
	 Vec vectildeqf = *dp;

    double scale = myPhase.pp->rho_r * 10 * myPhase.pp->l0;

    if (printpressure){

        // print converted pressure

        FILE * filesp = fopen(GetFilename(spname, mark), "w");
	     FILE * filedp = fopen(GetFilename(dpname, mark), "w");

        for (int j=0; j<N_; j++){
        for (int i=0; i<M_; i++){

            int dof = j*M_ + i;

            double phi = average_poro.at(dof);
            double invephi = 0.0;

            double tildeqf = 0.0;
				double q = 0.0;

            PetscCall(VecGetValues(vecq, 1, &dof, &q));
            PetscCall(VecGetValues(vectildeqf, 1, &dof, &tildeqf));

            vertex ccenter = cellcenter.at(dof);

            double add = myPhase.pp->rho_f * 10 * ccenter[1] * myPhase.pp->l0;

            if (phi > 1e-14){
                invephi = 1.0/sqrt(phi);
				}

            double qf = tildeqf * invephi;
				double qs = qf - 1.0/(1-phi) * (qf-q);

            fprintf(filesp, "%e ", qs);
            fprintf(filedp, "%e ", qf);
		      }
		  }

        fclose(filesp);
        fclose(filedp);

	 } else {

        // Print untreated pressure potentials
        printCellScalar(sp, spname, mark);
        printCellScalar(dp, dpname, mark);

	 }
   
    return 1;
}

int couple::printdivmass(int mark, const char * fieldname, 
                                   const vector<vertex>& stokesvel,
                                   const vector<vertex>& darcyvel){

    double massconv = 0.0;

    FILE * file = fopen(GetFilename(fieldname, mark), "w");

    int tolvert = N_*(M_+1)*3;

    vector<vertex> direction {{-1,0}, {0,-1}, {1,0}, {0,1}};

    // Phase averaged mass conservation velocity

    const valarray<double>& gwe = GaussWeightsEdge;
    const valarray<double>& gpe = GaussPointsEdge;

	 // cheat a bit
    double len = L_/(double)M_;

cout << "print mass conservation." << endl;
    for (int j=0; j<N_; j++){
    for (int i=0; i<M_; i++){
        massconv = 0.0;

        int left  = (j*(M_+1)+i)*3;
        int right = (j*(M_+1)+i+1)*3;

        int bottom = tolvert + (j*M_ + i)*3;
        int top    = tolvert + ((j+1)*M_ + i)*3;

        vector<int> dof = {left, bottom, right, top};

        for (int e=0; e<4; e++){
        for (int g=0; g<3; g++){

            vertex vel = darcyvel.at(dof[e]+g) *  edgeporo.at(dof[e]+g) + stokesvel.at(dof[e]+g);

            massconv += len/2.0*(vel[0]*direction[e][0] + vel[1]*direction[e][1]) *gwe[g];

		  }}
        if (abs(massconv) > 10e-14){
            // Record the cells not mass conserved
            cout << "The cell ("<< i << ", " << j << "), is not conserved." << endl;
				printf("The value of the integral of divergence is %.15f .\n", massconv);
		  }

        fprintf(file, "%.15f ", massconv);
    }}

    fclose(file);

    return 1;
}

/*
int couple::printsepmass(int mark, const char * name1, 
					                    const char * name2, 
											  const vector<vertex>& stokesvel,
											  const vector<vertex>& darcyvel,
											  const vector<double>& stokesp,
											  const vector<double>& darcyp){

    double mass1 = 0.0;
    double mass2 = 0.0;

    FILE * file = fopen(GetFilename(name1, mark), "w");
    FILE * file = fopen(GetFilename(name2, mark), "w");

    int tolvert = N_*(M_+1)*3;

    vector<vertex> direction {{-1,0}, {0,-1}, {1,0}, {0,1}};

    // Phase averaged mass conservation velocity

    const valarray<double>& gwe = GaussWeightsEdge;
    const valarray<double>& gpe = GaussPointsEdge;

    const valarray<double>& gwf = GaussWeightsFace;
    const vector<vertex>&   gpf = GaussPointsFace;

	 // cheat a bit
    double len = L_/(double)M_;

    double temp1 = 0.0;

    for (int j=0; j<N_; j++){
    for (int i=0; i<M_; i++){

        mass1 = 0.0;

        int left  = (j*(M_+1)+i)*3;
        int right = (j*(M_+1)+i+1)*3;

        int bottom = tolvert + (j*M_ + i)*3;
        int top    = tolvert + ((j+1)*M_ + i)*3;

        int celldof = j*M_ + i;

        vector<int> dof = {left, bottom, right, top};

        // Integration over cell of velocity
        for (int e=0; e<4; e++){
        for (int g=0; g<3; g++){

            vertex vel = stokesvel.at(dof[e]+g);

            massconv += len/2.0*(vel[0]*direction[e][0] + vel[1]*direction[e][1]) *gwe[g];

		  }}

        // Integration of pressure 
        for (int gg=0; gg<(int)gpf.size(); gg++){
            double pointporo = cellporo.at(celldof*gpf.size() + gg);
            temp1 += ;
		  }

    }}

    fclose(file);

    return 1;
}

*/
int couple::examineFullPorosity(int t){

    cout << "For time step " << t << endl;

    int tolvert = N_*(M_+1)*3;

    for (int j=0; j<N_; j++){
    for (int i=0; i<N_; i++){

        cout << "For cell (" << i << ", " << j <<"): " << endl; 

        // check cell centered porosity
        printf("cell-averaged : %.12f \n", average_poro.at(j*M_+i)); 

        // check surrounding porosity on gauss points
        cout << "Bottom edge : " ;
	     for (int g=0; g<3; g++){
            int bottom = tolvert + (j*M_+i)*3;
            printf("%.12f ", edgeporo.at(bottom+g));
	     }cout << endl;

        cout << "Top edge : " ;
	     for (int g=0; g<3; g++){
            int top = tolvert + ((j+1)*M_+i)*3;
            printf("%.12f ", edgeporo.at(top+g));
	     }cout << endl;
 
        cout << "Left edge : " ;
 	     for (int g=0; g<3; g++){
            int left = (j*(M_+1)+i)*3;
				printf("%.12f ", edgeporo.at(left+g));
	     }cout << endl;

        cout << "Right edge : " ;
	     for (int g=0; g<3; g++){
            int right = (j*(M_+1)+i+1)*3;
				printf("%.12f ", edgeporo.at(right+g));
	     }cout << endl;

    cout << endl; 

	 }cout << endl;}

    cout << endl;

    return 1;
}
