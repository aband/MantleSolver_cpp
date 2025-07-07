#include "driver.h"
#include "print.h"

/**
 * Customized transport function
 * Transporting phi with customized melting rate directly
 */

inline double R(const double& phi){

    double work = (3+phi-4*phi*phi)/3.0*phi;
    work = 1.0/sqrt(work);

    return work;
}

const double trueSoln(const MeshInfo& mi, const vertex& point, const indice& global, PhysProperty * pp, const double& L){

    // Constant porosity
    double work = 0.0;

    double phi = AssignPorosity(point, pp);

    work = -1.0*phi*phi*(1-phi)*(1-cosh(R(phi)*point[1])/cosh(R(phi)*L));

    return work;
}

const double trueSoln_p(const MeshInfo& mi, const vertex& point, const indice& global, PhysProperty * pp, const double& L){

    // Pievewise Constant porosity
    double work = 0.0;

    double phi = AssignPorosity(point, pp);

    double z = point[1];

    if (z<0){
        double r= R(phi);
        work = -1*phi*phi*(1-phi);

        double a = -1;
		  double b = (1-cosh(r*L))/sinh(r*L);

        work *= (1+a*cosh(r*z) + b*sinh(r*z));
    }

    return work;
}

inline vector<double> infcoeff(const int& n, const double& phi){

    vector<double> work;

    if (n<3){

	 work.resize(3);
    } else {
        work.resize(n);
    }

    work.at(0) = 1.0/(4*phi - 1);  
    work.at(1) = (-1.0*phi - 4*phi*phi*work.at(0))/(18*phi - 1);
    work.at(2) = (80.0/3.0 * phi*phi*phi*work.at(0) - 10*phi*phi*work.at(1))/(40*phi - 1);

    int nn = 0;

    for (int i=3; i<n; i++){
        nn = i*2;

        double c0 = 4.0/3.0 * phi*phi*phi * ((nn-4)*(nn+1)+ 4 + 4*nn);
        double c1 = 1.0/3.0 * phi*phi * ((nn-2)*(nn+3) + 4 + 2*(nn+2));
        double c2 = phi*((nn+5)*nn + 4) - 1.0;

        work.at(i) = (c0*work.at(i-2) - c1*work.at(i-1))/c2;
    }

    return work;
}

inline vector<double> recurcoeff(const int& n, const double& phi, 
					                  const int& r, const int& add){

    vector<double> work;

    assert(n>1);

    work.resize(n);

    vector<double> tmp;
    tmp.resize(n+1);

    // c1
    work.at(0) = 1.0;

    tmp.at(0) = 0.0;
    tmp.at(1) = 1.0;

    int nn = 0;

    for (int i=1; i<n; i++){
        nn = i*2+add;

        double c0 = 4.0/3.0 * phi*phi*phi * ((nn+r-4)*(nn+r+1)+ 4 + 4*(nn+r));
        double c1 = 1.0/3.0 * phi*phi * ((nn+r-2)*(nn+r+3) + 4 + 2*(nn+r+2));
        double c2 = phi*((nn+r+5)*(nn+r) + 4) - 1.0;

        work.at(i) = (c0*tmp.at(i-2+1) - c1*tmp.at(i-1+1))/c2;
        tmp.at(i+1) = work.at(i);

    }

    return work;
}

inline vector<double> newrecurcoeff(const int& n, const double& phi,
                                    const int& r, const int& add){

    vector<double> work;

    assert(n>0);

    work.resize(n);

    vector<double> tmp;

    tmp.resize(n+1);

    work.at(0) = 1.0;

    tmp.at(0) = 0.0;
	 tmp.at(1) = 1.0;

    int nn = 0;

    for (int i=1; i<n; i++){
        nn = i*2+add;

        double cn4 = 4.0/3.0 * phi*phi*phi*(nn+r-4)*(nn+r-3);
        double cn2 = 1.0/3.0 * phi*phi*(nn+r-2)*(nn+r-3);
        double cn  = phi*(nn+r)*(nn+r-3) - 1.0;

        work.at(i) = (cn4*tmp.at(i-2+1) - cn2*tmp.at(i-1+1))/cn;
        tmp.at(i+1)= work.at(i);

    }

    return work;
}

inline vector<double> newrecurcoeff_p(const int& n, const double& phi,
                                      const bool& nocut){

    vector<double> work;
    if (n<3){
        work.resize(3);
	 } else {
        work.resize(n);
    }

    // c4
    work.at(0) = phi*phi/(4*phi - 1);
    // c6
    work.at(1) = -1*(phi*phi*phi + 4*phi*phi*work.at(0)) / (18*phi - 1);
	 // c8
    work.at(2) = (80.0/3.0*phi*phi*phi*work.at(0) - 10*phi*phi*work.at(1))/
				(40*phi-1);

    int nn = 0;
    if (n>3){

        if (nocut){
            for (int i=3; i<n; i++){

                // start with c10
                nn = i*2 + 4; 

                double cn4 = 4.0/3.0 * phi*phi*phi*(nn-4)*(nn-3);
                double cn2 = 1.0/3.0 * phi*phi*(nn-2)*(nn-3);
                double cn  = phi*(nn)*(nn-3) - 1.0;

                work.at(i) = (cn4*work.at(i-2) - cn2*work.at(i-1))/cn;
 
            }
        } else {
            for (int i=3; i<n; i++){

                work.at(i) = 0.0;
            }
	     }
    }

    return work;
}

inline double recurcoeff_derive(const vector<double>& infcoeff,
                         const vector<double>& ecoeff,
		                   const int& r, const int& cutoff,
								 const double& phi, const double& z, const double& c){

    double work = 0.0;

    for (int i=0; i<cutoff; i++){
        work += infcoeff.at(i)*phi*phi*(4+i*2)*pow(abs(z),4+i*2-1)
                + c* ecoeff.at(i)*phi*phi*(4+i*2+r)*pow(abs(z),4+i*2+r-1);
    }
    return work;
}

inline double recurcoeff_inte(const vector<double>& infcoeff,
                              const vector<double>& ecoeff,
		                        const int& r, const int& cutoff,
								      const double& phi, const double& z, 
										const double& c){

    double work = 0.0;

    for (int i=0; i<cutoff; i++){
        work += infcoeff.at(i)*pow(abs(z),i*2+1)/(i*2+1.0)
             + c*ecoeff.at(i)*pow(abs(z),i*2+r+1)/(i*2+1.0+r);
    }

    return work;
}

/*
const double trueSoln_q(const MeshInfo& mi, const vertex& point, const indice& global, PhysProperty * pp, const double& L){

    double work = 0.0;
    double work1 = 0.0;

    double phi = AssignPorosity(point, pp);
    phi = 0.01;

    double z = point[1];

    double r1 = (3+sqrt(9+4/phi))/2;
    double r2 = (-5+sqrt(9+4/phi))/2;

    int cutoff = 5;
    vector<double> coeff = infcoeff(cutoff,phi);
    vector<double> ocoeff = recurcoeff(cutoff,phi,r2,1);
    vector<double> ecoeff = recurcoeff(cutoff,phi,r2,0);

    vector<double> c_p = newrecurcoeff_p(cutoff, phi,true); 
    vector<double> odd_c_h = newrecurcoeff(cutoff, phi, r1, 1);
    vector<double> even_c_h = newrecurcoeff(cutoff, phi, r1, 0);

    double p1 = 0.0; 
    double a1 = 0.0;
    double b1 = 0.0;

    double newp1 = 0.0;
    double newa1 = 0.0;
    double newb1 = 0.0;

    for (int i=0; i<cutoff; i++){
        p1 += coeff.at(i)*pow(2,i*2);
        a1 += ecoeff.at(i)*pow(2,i*2+r2);
        b1 += ocoeff.at(i)*pow(2,i*2+1+r2);

        newp1 += c_p.at(i)*pow(2,i*2+4); 
        newa1 += even_c_h.at(i)*pow(2,i*2+r1);
        newb1 += odd_c_h.at(i)*pow(2,i*2+1+r1); 
	//	  cout << c_p.at(i) << endl;
//		  cout << coeff.at(i) << endl;
    }

    double p2 = 0.0; 
    double a2 = 0.0;
    double b2 = 0.0;

    for (int i=0; i<cutoff; i++){
        p2 += coeff.at(i)*pow(1.2,i*2+4);
        a2 += ecoeff.at(i)*pow(1.2,i*2+r2);
        b2 += ocoeff.at(i)*pow(1.2,i*2+1+r2);
    }

    double c1 = 0.0;
    double c2 = 0.0;

    double newc1 = 0.0;
    double newc2 = 0.0;

    double multi = 1.0/(a1*b2-a2*b1);

    c1 = b1*p2-b2*p1;
    c2 = a2*p1-a1*p2;

    c1*=multi;
    c2*=multi;

    c1 = -p1/a1;
    newc1 = -newp1/newa1;

    if (z<0){
        double scale = phi*phi*z*z*z*z;


//        work = scale * 1.0/(4*phi-1)  
//                + phi*phi/(1-4*phi)*pow(L,4-r1)*pow(-1*z,r1);


        for (int i=0; i<cutoff; i++){
            work1 += coeff.at(i)*pow(z,i*2) * scale;
				//		  + c1*ecoeff.at(i)*pow(abs(z),i*2+r2);
                 // + c2*ocoeff.at(i)*pow(abs(z),i*2+1+r2));
            work += c_p.at(i)*pow(z,i*2+4)
                  + newc1*even_c_h.at(i)*pow(abs(z),i*2+r1);
//cout << work << "  " << work1 << endl;
        }
        //work *= scale;

    }

    return work;
}
*/

const double trueSoln_q(const MeshInfo& mi, const vertex& point, const indice& global, PhysProperty * pp, const double& L){

    double work = 0.0;

    double phi = 0.01;

    double z = point[1];
    double r = (3+sqrt(9+4/phi))/2;
    double rn = (3-sqrt(9+4/phi))/2;

    int cutoff = 5;
    vector<double> c_p = newrecurcoeff_p(cutoff, phi, true); 
    vector<double> odd_c_h = newrecurcoeff(cutoff, phi, r, 1);
    vector<double> even_c_h = newrecurcoeff(cutoff, phi, r, 0);
    vector<double> even_c_hn = newrecurcoeff(cutoff, phi, rn, 0);

    double p1 = 0.0; 
    double a1 = 0.0;
    double b1 = 0.0;
    double p2 = 0.0;
    double a2 = 0.0;
    double b2 = 0.0;
    double d1 = 0.0;
    double d2 = 0.0;

    for (int i=0; i<cutoff; i++){
        p1 += c_p.at(i)*pow(0.2,i*2+4); 
        a1 += even_c_h.at(i)*pow(0.2,i*2+r); 
        d1 += even_c_hn.at(i)*pow(0.2,i*2+rn);
        //b1 += odd_c_h.at(i)*pow(0.2,i*2+1+r); 

        p2 += c_p.at(i)*pow(2,i*2+4); 
        a2 += even_c_h.at(i)*pow(2,i*2+r); 
        d2 += even_c_hn.at(i)*pow(2,i*2+rn);
        //b2 += odd_c_h.at(i)*pow(2,i*2+1+r); 

        //p2 += c_p.at(i)*pow(0.0001,i*2+4);
        //a2 += even_c_h.at(i)*pow(0.0001,i*2+r);
        //b2 += odd_c_h.at(i)*pow(0.0001,i*2+1+r);
    }	

    //double c2 = -p1/b1;
    double multi = 1.0/(a1*d2-a2*d1);
    double c1 = d1*p2-d2*p1;
    double c2 = a2*p1-a1*p2;

    c1*=multi;
    c2*=multi;

    c1 = -p2/a2;
    c2 = 0.0;

    if (z<-0.2){
        for (int i=0; i<cutoff; i++){
             work += c_p.at(i)*pow(abs(z),i*2+4)
                  + c1*even_c_h.at(i)*pow(abs(z),i*2+r)
                  + c2*even_c_hn.at(i)*pow(abs(z),i*2+rn);
                  //+ c2*odd_c_h.at(i)*pow(abs(z),i*2+r+1);
        }
    }

    return work;
}
// True solution for pressure potentials

const std::array<double,2> trueSolnq(const MeshInfo& mi, const vertex& point, const indice& global, PhysProperty * pp, const double& L){

    std::array<double,2> work {0.0,0.0};

    double phi = AssignPorosity(point, pp);
    double r = R(phi);
    double z = point[1];

    work.at(0) = (1 - phi) * (z - 1.0/r * (sinh(r*z)/cosh(r*L)));
    work.at(1) = (1 - phi) * (z - (1-4*phi)/(3+phi-4*phi*phi) * phi/r * (sinh(r*z)/cosh(r*L)));

    return work;
}

const std::array<double,2> trueSolnq_p(const MeshInfo& mi, const vertex& point,
const indice& global, PhysProperty * pp, const double& L){

    std::array<double,2> work {0.0,0.0};

    double phi = AssignPorosity(point, pp);
    double r = R(phi);
    double z = point[1];

    double b= 0.0;
    double c= 0.0;

    if (z<0){
        double r= R(phi);
        double tmp = -1*phi*phi*(1-phi);

        double a = -1;
		  b = (1-cosh(r*L))/sinh(r*L);
        c = -b * (1-phi)/r;

        tmp *= (1+a*cosh(r*z) + b*sinh(r*z));
        work[1] = (1-phi)*z;


        work.at(0) = (1 - phi)*(z - 1/r * sinh(r*z) +b/r * cosh(r*z));

        //work.at(1) = (1 - phi) * (z - (1-4*phi)/(3+phi-4*phi*phi) * phi/r * (sinh(r*z)/cosh(r*2)));

        work.at(1) = work.at(0) - phi*(1-phi)*(-1*r*sinh(r*z) + b*r*cosh(r*z));
    }else {
        work[1] = z ;
		  work[0] = 0.0;
	 } 

    return work;
}

inline double newrecurcoeff_deriv(const vector<double>& c_p,
                                  const vector<double>& even_c_h, 
                                  const vector<double>& even_c_hn,
											 const double& c1,
											 const double& c2,
											 const int& r1,
											 const int& r2,
											 const double& phi, 
											 const double& z,
											 const int& cutoff){

    double work = 0.0;

    for (int i=0; i<cutoff; i++){
        work += c_p.at(i)*pow(abs(z),2*i+4-1)*(2*i+4)+
					 c1*even_c_h.at(i)*pow(abs(z),i*2+r1-1)*(2*i+r1)+
					 c2*even_c_hn.at(i)*pow(abs(z),i*2+r2-1)*(2*i+r2);
    }

    return work;
}

inline double newrecurcoeff_inte(const vector<double>& c_p,
                              const vector<double>& even_c_h,
										const vector<double>& even_c_hn,
										const double& c1,
										const double& c2,
									   const int& r1,
										const int& r2,
										const double& phi, 
										const double& z,
										const int& cutoff){

    double work = 0.0;

	 for (int i=0; i<cutoff; i++){
        work += c_p.at(i)*pow(abs(z),2*i+1)/(2*i+1)/phi/phi+
		   		 c1*even_c_h.at(i)*pow(abs(z),i*2+r1-4+1)/(2*i+r1-4+1)/phi/phi+
	  	   		 c2*even_c_hn.at(i)*pow(abs(z),i*2+r2-4+1)/(2*i+r2-4+1)/phi/phi;
    } 

    return work;
}

const std::array<double,2> trueSolnq_q(const MeshInfo& mi, const vertex& point, const indice& global, PhysProperty * pp, const double& L){

    std::array<double,2> work {0.0,0.0};

    double z = point[1];
    double phi = 0.01;
    double r1 = (3+sqrt(9+4/phi))/2;
    double r2 = (3-sqrt(9+4/phi))/2;

    int cutoff = 5;
    vector<double> c_p = newrecurcoeff_p(cutoff, phi, true); 
 
    vector<double> even_c_h = newrecurcoeff(cutoff, phi, r1, 0);
    vector<double> even_c_hn = newrecurcoeff(cutoff, phi, r2, 0);

    double p1 = 0.0; 
    double a1 = 0.0;
    double p2 = 0.0;
    double a2 = 0.0;
    double d1 = 0.0;
    double d2 = 0.0;

    for (int i=0; i<cutoff; i++){
        p1 += c_p.at(i)*pow(0.2,i*2+4); 
        a1 += even_c_h.at(i)*pow(0.2,i*2+r1); 
        d1 += even_c_hn.at(i)*pow(0.2,i*2+r2);

        p2 += c_p.at(i)*pow(2,i*2+4); 
        a2 += even_c_h.at(i)*pow(2,i*2+r1); 
        d2 += even_c_hn.at(i)*pow(2,i*2+r2);

        //p2 += c_p.at(i)*pow(0.0001,i*2+4);
        //a2 += even_c_h.at(i)*pow(0.0001,i*2+r);
        //b2 += odd_c_h.at(i)*pow(0.0001,i*2+1+r);
    }	

    double multi = 1.0/(a1*d2-a2*d1);
    double c1 = d1*p2-d2*p1;
    double c2 = a2*p1-a1*p2;

    c1*=multi;
    c2*=multi;

    c1=-p2/a2;
    c2= 0.0;

    double tmp2 = newrecurcoeff_deriv(c_p, even_c_h, even_c_hn, c1,c2,
                              r1, r2, phi, z, cutoff); 
    double tmp1 = newrecurcoeff_inte(c_p, even_c_h, even_c_hn, c1,c2,
                              r1, r2, phi, z, cutoff); 
    if (z < -0.0){

        //work.at(0) = tmp1;

        work.at(1) = z - 1.0/3.0*phi*pow(z,3) + tmp2*(1-4*phi*z*z)/3.0;

        //work.at(0) = work.at(1) + tmp2/phi/z/z;
        work.at(0) = tmp1;
    }else{
        work.at(1) = z;
        work.at(0) = 0.0;
    }

    return work;
} 

/*
const std::array<double,2> trueSolnq_q(const MeshInfo& mi, const vertex& point, const indice& global, PhysProperty * pp, const double& L){

    std::array<double,2> work {0.0,0.0};

    double phi = 0.01;

    double z = point[1];

    double r1 = (3+sqrt(9+4/phi))/2;
    double r2 = (3-sqrt(9+4/phi))/2;

    double tmp1 = phi*phi/(1-4*phi) * (r1*pow(L,4-r1)*pow(abs(z),r1-1) - 4*pow(abs(z),3));

    int cutoff = 10;
    double r = (-5+sqrt(9+4/phi))/2;

    vector<double> coeff = infcoeff(cutoff,phi);
    vector<double> ecoeff = recurcoeff(cutoff,phi,r,0);

    double p1 = 0.0;
	 double a1 = 0.0;

    for (int i=0; i<cutoff; i++){
        p1 += coeff.at(i)*pow(2,i*2);
        a1 += ecoeff.at(i)*pow(2,i*2+r);
    }

    double c1 = -p1/a1;

    //cout << c1 << "  " << pow(L,4-r1)/(1-4*phi) << "  " << 1.0/ pow(2,r1-4);

    double tmp2 = recurcoeff_derive(coeff, ecoeff, r, cutoff, phi, z, c1); 

    // ql
	 double tt = 0.0;
	 double ttt = 0.0; 

    int N = mi.MPIglobalCellSize[1];
    if (z < 0.2){
        tt = 1.0/(1-4*phi) * (z + pow(L,4-r1)*pow(abs(z),r1-3)/(r1-3));

        //work.at(1) = tt - tmp1/phi/z/z;

        ttt = z - 1.0/3.0 * phi * pow(z,3) - tmp1 * (1-4*phi*z*z)/3.0;

        work.at(0) = recurcoeff_inte(coeff,ecoeff,r,cutoff,phi,z,c1);
		  //+ c1*(phi*phi*(N-10)/60.0*pow(abs(z),r+11)); 
		  //+ 1.3*c1*(phi*(N+20)/80*pow(abs(z),r+2)); 

        work.at(1) = z - 1.0/3.0 * phi * pow(z,3) + tmp2 * (1-4*phi*z*z)/3.0; 
		  //+ c1*(phi*phi*0.04*(N-20)/60.0*pow(abs(z),r+6)); 
		  //+ c1*(phi*phi*0.1*(N)/60.0*pow(abs(z),r+6)); 

    } else {

        work.at(1) = point[1];

        work.at(0) = 0.0;
    }
    return work;
}
*/

int printExactPorosity(int mark, const MeshInfo& mi, const char * fieldname, PhysProperty * pp){

    char * filename = (char *)malloc(strlen(fieldname)+10+4);

    char n_char[10];
    std::sprintf(n_char,"%d",mark);
    strcpy(filename, fieldname);
    strcat(filename, n_char);
    strcat(filename, ".dat");

    FILE * sol = fopen(filename,"w");

//    FILE * exactql = fopen("exactql", "w");
//    FILE * exactqs = fopen("exactqs", "w");

    for(int j=0; j<mi.MPIglobalCellSize[1]; j++){
    for(int i=0; i<mi.MPIglobalCellSize[0]; i++){

        vertex local {0.0, 0.0};

        vector<vertex> corners = extractCorners(mi, {i,j});

        vertex global = GaussMapPointsFace(local, corners);

        fprintf(sol, "%.16f ", AssignPorosity(global,pp));

//        std::array<double,2> exactvals = trueSolnq(mi, global, {i,j}, pp, 2);

//        fprintf(exactql, "%.16f ", -1.0*exactvals[0]);
//        fprintf(exactqs, "%.16f ", -1.0*exactvals[1]);

    }fprintf(sol, "\n");}

    fclose(sol);
//    fclose(exactql);
//    fclose(exactqs);

    return 1;
}

/**
 * Print velocity on edge gauss points
 */
int Driver::printVelEdgeGauss_case(int mark, PhysProperty * pp){

    FILE * dvx = fopen(GetFilename("darcyvelx", mark),"w");
    FILE * dvy = fopen(GetFilename("darcyvely", mark),"w");

    FILE * svx = fopen(GetFilename("stokesvelx", mark),"w");
    FILE * svy = fopen(GetFilename("stokesvely", mark),"w");

    FILE * phasevx = fopen(GetFilename("phasevx",mark),"w");
    FILE * phasevy = fopen(GetFilename("phasevy",mark),"w");

    FILE * gaussgridx = fopen("gaussgridx.dat", "w");
    FILE * gaussgridy = fopen("gaussgridy.dat", "w");

    FILE * exactv = fopen("exactv.dat", "w");

    const valarray<double>& gwe = GaussWeightsEdge;
    const valarray<double>& gpe = GaussPointsEdge;

    std::vector<vertex> gaussp;
    gaussp.resize(gpe.size());

    vertexSet edge;

    vector<vertex> darcyvel;  darcyvel.resize(gaussp.size());
    vector<vertex> stokesvel; stokesvel.resize(gaussp.size());

    for (int j=0; j<mi.MPIglobalCellSize[1]; j++){
    for (int i=0; i<mi.MPIglobalCellSize[0]; i++){

        indice gcell {i,j};
        indice gcellout;
        vertexSet corners = extractCorners(mi, gcell);

        darcyvel.clear(); darcyvel.resize(gaussp.size());
        stokesvel.clear(); stokesvel.resize(gaussp.size());
        
        // Horizontal edges only
        edge = {corners.at(0), corners.at(1)};

        for (int g=0; g<gpe.size(); g++){

            gaussp.at(g) = GaussMapPointsEdge({gpe[g]}, edge);

        }


        vector<vertex> vel_relative = 
        ExtractVelocity(&sresult_->vel_darcy, &sresult_->g_darcy,
                    refArrayDarcyEssen_,mi,
                    gaussp, gcell,*hdiv_,*basis_,{1});
    
        vector<vertex> vel_stokes = 
        ExtractVelocity(&sresult_->vel_stokes, &sresult_->g_stokes,
                    refArrayStokesEssen_,mi,
                    gaussp, gcell,*br_,*basis_,{1});

        for (int g=0; g<gpe.size(); g++){

//            gaussp.at(g) = GaussMapPointsEdge({gpe[g]}, edge);
/*
        vector<vertex> vel_relative = 
        ExtractVelocity(&sresult_->vel_darcy, &sresult_->g_darcy,
                    refArrayDarcyEssen_,mi,
                    gaussp, gcell,*hdiv_,*basis_,{1});
    
        vector<vertex> vel_stokes = 
        ExtractVelocity(&sresult_->vel_stokes, &sresult_->g_stokes,
                    refArrayStokesEssen_,mi,
                    gaussp, gcell,*br_,*basis_,{1});
*/

            darcyvel.at(g) = AssignPorosity(gaussp.at(g),pp) * vel_relative.at(g);

            stokesvel.at(g) = vel_stokes.at(g);

            fprintf(dvx, "%e ", darcyvel.at(g)[0]);
            fprintf(dvy, "%e ", darcyvel.at(g)[1]);

            fprintf(svx, "%e ", stokesvel.at(g)[0]);
            fprintf(svy, "%e ", stokesvel.at(g)[1]);

            fprintf(gaussgridx, "%e ", gaussp.at(g)[0]);
            fprintf(gaussgridy, "%e ", gaussp.at(g)[1]);

            fprintf(exactv, "%e ", trueSoln_q(mi, gaussp.at(g), gcell, pp, mi.H/2.0));

            fprintf(phasevx, "%e ", stokesvel.at(g)[0] + darcyvel.at(g)[0]);
            fprintf(phasevy, "%e ", stokesvel.at(g)[1] + darcyvel.at(g)[1]);
        }

     }}

     // Add top boundary
     for (int i=0; i<mi.MPIglobalCellSize[0]; i++){

        int j= mi.MPIglobalCellSize[1]-1;

        indice gcell {i,j};

        vertexSet corners = extractCorners(mi, gcell);

        edge = {corners.at(2), corners.at(3)};

        for (int g=0; g<gpe.size(); g++){

            gaussp.at(g) = GaussMapPointsEdge({gpe[g]}, edge);

        vector<vertex> vel_relative = 
        ExtractVelocity(&sresult_->vel_darcy, &sresult_->g_darcy,
                    refArrayDarcyEssen_,mi,
                    gaussp, gcell,*hdiv_,*basis_,{1});
    
        vector<vertex> vel_stokes = 
        ExtractVelocity(&sresult_->vel_stokes, &sresult_->g_stokes,
                    refArrayStokesEssen_,mi,
                    gaussp, gcell,*br_,*basis_,{1});

            darcyvel.at(g) = AssignPorosity(gaussp.at(g),pp) * vel_relative.at(g);
            stokesvel.at(g) = vel_stokes.at(g);

            fprintf(dvx, "%e ", darcyvel.at(g)[0]);
            fprintf(dvy, "%e ", darcyvel.at(g)[1]);

            fprintf(svx, "%e ", stokesvel.at(g)[0]);
            fprintf(svy, "%e ", stokesvel.at(g)[1]);

            fprintf(gaussgridx, "%e ", gaussp.at(g)[0]);
            fprintf(gaussgridy, "%e ", gaussp.at(g)[1]);

            fprintf(exactv, "%e ", 0.0);

            fprintf(phasevx, "%e ", stokesvel.at(g)[0] + darcyvel.at(g)[0]);
            fprintf(phasevy, "%e ", stokesvel.at(g)[1] + darcyvel.at(g)[1]);
           
        }

     }

/*
     }fprintf(dvx, "\n ");
      fprintf(dvy, "\n ");
      fprintf(svx, "\n ");
      fprintf(svy, "\n ");
      fprintf(gaussgridx, "\n ");
      fprintf(gaussgridy, "\n ");}
*/
      fclose(dvx);
      fclose(dvy);
      fclose(svx);
      fclose(svy);
      fclose(gaussgridx);
      fclose(gaussgridy);
      fclose(exactv);
      fclose(phasevx);
		fclose(phasevy);

    return 1;
}

// Lp norm for velocity error
double Driver::errorNorm(int mark, PhysProperty * pp, int norm){

    // Compute error norm 
    double totalerror = 0.0;

    const valarray<double>& gwf = GaussWeightsFace;
    const vector<vertex>& gpf = GaussPointsFace;

    std::vector<vertex> gaussp;
    gaussp.resize(gpf.size());

    vector<vertex> darcyvel;  darcyvel.resize(gaussp.size());
    vector<vertex> stokesvel; stokesvel.resize(gaussp.size());

    for (int j=0; j<mi.MPIglobalCellSize[1]; j++){
    for (int i=0; i<mi.MPIglobalCellSize[0]; i++){
        indice gcell {i,j};

        vertexSet corners = extractCorners(mi, gcell);

        for (unsigned int g=0; g<gwf.size(); g++){
            gaussp.at(g) = GaussMapPointsFace(gpf[g],corners);
        }

        vector<vertex> vel_relative = 
        ExtractVelocity(&sresult_->vel_darcy, &sresult_->g_darcy,
                    refArrayDarcyEssen_,mi,
                    gaussp, gcell,*hdiv_,*basis_,{1});
    
        vector<vertex> vel_stokes = 
        ExtractVelocity(&sresult_->vel_stokes, &sresult_->g_stokes,
                    refArrayStokesEssen_,mi,
                    gaussp, gcell,*br_,*basis_,{1});

        double work = 0.0;
        for (unsigned int g=0; g<gwf.size(); g++){
            double jac = abs(GaussJacobian(gpf[g],corners));
            double gw = gwf[g];

            darcyvel.at(g) = AssignPorosity(gaussp.at(g),pp) * vel_relative.at(g);
            stokesvel.at(g) = vel_stokes.at(g);

            work += jac * gw * pow(abs(stokesvel.at(g)[1])-abs(trueSoln_q(mi,gaussp.at(g),gcell,pp,mi.H/2.0)),norm);

        }

        double area = mi.cellArea.at(FlatIndic(mi,gcell));

        totalerror += work;
    }}

    return sqrt(totalerror);
}

// Lp norm for pressure error
/*
std::array<double,2> Driver::errorP(int mark, PhysProperty * pp, int norm){

    double add1 = 0.0, add2 = 0.0;

	 std::array<double, 2> work {0.0,0.0};

    std::array<double, 2> truesol {0.0,0.0};

    for (int j=0; j<mi.MPIglobalCellSize[1]; j++){
    for (int i=0; i<mi.MPIglobalCellSize[0]; i++){
        indice gcell {i,j};

        vertexSet corners = extractCorners(mi, gcell);

        // Get middle point
        vertex local {0.0,0.0};	
        vertex mid = GaussMapPointsFace(local, corners);

        double area = mi.cellArea.at(FlatIndic(mi,gcell));

        // Get true solution
        truesol = trueSolnq(mi, mid, {i,j}, pp, 2);


        work.at(0) += add1;
        work.at(1) += add2;
    }}

    return work;
}
*/

double averagePhi(PhysProperty * pp, int i, int j, const MeshInfo& mi){

    const valarray<double>& gwf = GaussWeightsFace;
    const vector<vertex>&   gpf = GaussPointsFace;

    double work = 0.0;
    double area = 0.0;

    for (unsigned int g=0; g<gwf.size(); g++){
        vector<vertex> corners = extractCorners(mi, {i,j});
        vertex mapped = GaussMapPointsFace(gpf[g],corners);

        // Get HD and CD from reconstruction at this gaussian point
        double phif = 0.0;

        phif = AssignPorosity(mapped, pp); 
//cout << phif << endl;
        // ======================================================

        double jac = abs(GaussJacobian(gpf[g],corners));
        double gw = gwf[g];
        work += gw * jac * phif;
        area += gw * jac; 
    }

    work /= area;

    return work;
}


int Driver::printSimplePressure_case(int mark, PhysProperty * pp){

    Vec vectildeqf;
    Vec vecq;

    PetscCall(VecNestGetSubVec(Result_->y, 0, &vecq));   
    PetscCall(VecNestGetSubVec(Result_->y, 1, &vectildeqf));

    FILE * fqs = fopen(GetFilename("qs", mark), "w");
    FILE * fqf = fopen(GetFilename("qf", mark), "w");

    FILE * fstokes = fopen(GetFilename("rawstokesq", mark), "w"); 
    FILE * fdarcy  = fopen(GetFilename("rawdarcyq", mark), "w");

    for (int j=0; j<mi.MPIglobalCellSize[1]; j++){
    for (int i=0; i<mi.MPIglobalCellSize[0]; i++){

        indice globalcell {i,j};
        int nelem = FlatIndic(mi, globalcell);
        double q, tildeqf;

        PetscCall(VecGetValues(vectildeqf, 1, &nelem, &tildeqf));
        PetscCall(VecGetValues(vecq, 1, &nelem, &q));

        // Get Porosity
        vertex local {0.0,0.0};

        basis_->GetCorners(mi, globalcell);

        vertex global = GaussMapPointsFace(local, basis_->corners());

        double phif = AssignPorosity(global, pp);
		  double avephi = averagePhi(pp, i, j, mi);
		  double coef = 0.0;
        // Adjust phif
        if (avephi > 1e-16) {
            coef = 1.0/sqrt(avephi);
        }

        // Reterive original physical variables with physical units
        double qf = tildeqf *coef;
        double qs = -qf - 1.0/(1-phif)*(-qf-q);

        fprintf(fqs, "%.16f ", qs);
        fprintf(fqf, "%.16f ", qf);
        fprintf(fstokes, "%.16f ", q);
        fprintf(fdarcy, "%.16f ", tildeqf);
    }fprintf(fqs, "\n");
     fprintf(fqf, "\n");
     fprintf(fstokes, "\n");
     fprintf(fdarcy, "\n");}

    return 1;
}

// Print averaged-shifted pressure potential
int Driver::printShiftedPressure_case(int mark, PhysProperty * pp){

    Vec vectildeqf;
    Vec vecq;

    PetscCall(VecNestGetSubVec(Result_->y, 0, &vecq));   
    PetscCall(VecNestGetSubVec(Result_->y, 1, &vectildeqf));

    FILE * fqs = fopen(GetFilename("qs", mark), "w");
    FILE * fqf = fopen(GetFilename("qf", mark), "w");

    FILE * exactql = fopen("exactql", "w");
    FILE * exactqs = fopen("exactqs", "w");

    FILE * errorqs_all = fopen("errorqs", "w");
    FILE * errorql_all = fopen("errorql", "w");

    double sumqs = 0.0;

    std::vector<double> qs_vec;
	 std::vector<double> ql_vec;

    std::vector<double> qs_exact;
    std::vector<double> ql_exact;

    std::vector<double> q_bar;
    std::vector<double> q_vec;

    std::vector<bool> mask;

    qs_vec.resize(mi.MPIglobalCellSize[0]*mi.MPIglobalCellSize[1]);
    ql_vec.resize(mi.MPIglobalCellSize[0]*mi.MPIglobalCellSize[1]);

    qs_exact.resize(mi.MPIglobalCellSize[0]*mi.MPIglobalCellSize[1]);
    ql_exact.resize(mi.MPIglobalCellSize[0]*mi.MPIglobalCellSize[1]);

    mask.resize(mi.MPIglobalCellSize[0]*mi.MPIglobalCellSize[1]);

    q_vec.resize(mi.MPIglobalCellSize[0]*mi.MPIglobalCellSize[1]);
double sumexact = 0.0;
double sumqs_plus = 0.0;
double sumqs_minus = 0.0;

//    int startj = mi.MPIglobalCellSize[1]/2;
    int startj = mi.MPIglobalCellSize[1];

    int halfsize = mi.MPIglobalCellSize[1]/2;

    //for (int j=0; j<mi.MPIglobalCellSize[1]; j++){
    for (int j=0; j<startj; j++){
    for (int i=0; i<mi.MPIglobalCellSize[0]; i++){

        indice globalcell {i,j};
        int nelem = FlatIndic(mi, globalcell);
        double q, tildeqf;

        PetscCall(VecGetValues(vectildeqf, 1, &nelem, &tildeqf));
        PetscCall(VecGetValues(vecq, 1, &nelem, &q));

        // Get Porosity
        vertex local {0.0,0.0};

        basis_->GetCorners(mi, globalcell);

        vertex global = GaussMapPointsFace(local, basis_->corners());

        // Pointwise porosity evaluated at the middle point
        double phif = AssignPorosity(global, pp);
      
        // Averaged porosity of the cell
        double avephi = averagePhi(pp, i, j, mi);
        double coef = 0.0;

        mask.at(nelem) = false;
        // Adjust phif
        if (avephi > 1e-16) {
            coef = 1.0/sqrt(avephi);
				mask.at(nelem) = true;

        }

        double ql = tildeqf *coef;
//        double qs = -ql + 1.0/(1-avephi)*(ql+q);
        double qs = -ql + 1.0/(1-phif)*(ql+q);

        ql_vec.at(nelem) = -ql;
        qs_vec.at(nelem) = qs;
        q_vec.at(nelem) = q;
//cout << mi.H/2.0 << endl;
        std::array<double,2> exactval = trueSolnq_q(mi, global, {i,j}, pp, mi.H/2.0);

        ql_exact.at(nelem) = -1.0*exactval[0];
        qs_exact.at(nelem) = -1.0*exactval[1];

        sumqs += qs;
sumexact += -1.0*exactval[1];
    }} 


    // ====================================================================
    // Separate with \phi=0 and \phi \neq 0
    double sum_exactqs_half1 = 0.0;
    double sum_exactqs_half2 = 0.0;

    double sum_exactql_half1 = 0.0;
    double sum_exactql_half2 = 0.0;

    double sum_vecqs_half1 = 0.0;
    double sum_vecqs_half2 = 0.0;

    double sum_vecql_half1 = 0.0;
    double sum_vecql_half2 = 0.0;

    for (int j=0; j<halfsize; j++){
    for (int i=0; i<mi.MPIglobalCellSize[0]; i++){

        sum_exactqs_half1 += qs_exact.at(FlatIndic(mi,{i,j}));
        sum_exactqs_half2 += qs_exact.at(FlatIndic(mi,{i,j+halfsize}));

        sum_exactql_half1 += ql_exact.at(FlatIndic(mi,{i,j}));
        sum_exactql_half2 += ql_exact.at(FlatIndic(mi,{i,j+halfsize}));

        sum_vecqs_half1 += qs_vec.at(FlatIndic(mi,{i,j}));
        sum_vecqs_half2 += qs_vec.at(FlatIndic(mi,{i,j+halfsize}));

        sum_vecql_half1 += ql_vec.at(FlatIndic(mi,{i,j}));
        sum_vecql_half2 += ql_vec.at(FlatIndic(mi,{i,j+halfsize}));
    }}

//    cout << halfsize << endl;

    sum_exactqs_half1 /= (double)mi.MPIglobalCellSize[0]*halfsize; 
    sum_exactqs_half2 /= (double)mi.MPIglobalCellSize[0]*halfsize; 

    sum_exactql_half1 /= (double)mi.MPIglobalCellSize[0]*halfsize; 
    sum_exactql_half2 /= (double)mi.MPIglobalCellSize[0]*halfsize; 

    sum_vecqs_half1 /= (double)mi.MPIglobalCellSize[0]*halfsize; 
    sum_vecqs_half2 /= (double)mi.MPIglobalCellSize[0]*halfsize; 

    sum_vecql_half1 /= (double)mi.MPIglobalCellSize[0]*halfsize; 
    sum_vecql_half2 /= (double)mi.MPIglobalCellSize[0]*halfsize; 

//    cout << sum_exactqs_half1 - sum_vecqs_half1 << "  " << sum_exactqs_half2 - sum_vecqs_half2 << endl;
//    cout << sum_vecqs_half1 << "  " << sum_vecqs_half2 << endl;

    double errorqs_half1 = 0.0;
    double errorql_half1 = 0.0;

    double errorqs_half2 = 0.0;
    double errorql_half2 = 0.0;

    double sum_errorqs = 0.0;
    double sum_errorql = 0.0;

    int N = halfsize*2;
//    int cutoff = halfsize*0.875;
    int cutoff = halfsize;

    // Renormalized error results
    for (int j=0; j<halfsize; j++){
    for (int i=0; i<mi.MPIglobalCellSize[0]; i++){

        double area1 = mi.cellArea.at(FlatIndic(mi,{i,j}));
        double area2 = mi.cellArea.at(FlatIndic(mi,{i,j+halfsize}));

//cout << diff_half1 << "  " << diff_half2 << endl;

        double diff_half1_ql = 0.0;
        double diff_half1    = 0.0;

        if (j > halfsize - cutoff){
        diff_half1_ql = ql_exact.at(FlatIndic(mi,{i,j})) - sum_exactql_half1 - 
                        (ql_vec.at(FlatIndic(mi,{i,j})) - sum_vecql_half1);

        diff_half1 = qs_exact.at(FlatIndic(mi,{i,j})) - sum_exactqs_half1 - 
                     (qs_vec.at(FlatIndic(mi,{i,j})) - sum_vecqs_half1);
        }

        double diff_half2_ql = 0.0;
        double diff_half2 = 0.0;

        diff_half2_ql = ql_exact.at(FlatIndic(mi,{i,j+halfsize})) - sum_exactql_half2 - 
                              (ql_vec.at(FlatIndic(mi,{i,j+halfsize})) - sum_vecql_half2);

        diff_half2 = qs_exact.at(FlatIndic(mi,{i,j+halfsize})) - sum_exactqs_half2 - 
                           (qs_vec.at(FlatIndic(mi,{i,j+halfsize})) - sum_vecqs_half2);

        errorqs_half1 += pow(diff_half1,2)*area1; 
        errorqs_half2 += pow(diff_half2,2)*area2; 

        sum_errorqs += pow(diff_half1,2)*area1 + pow(diff_half2,2)*area2;

        sum_errorql += pow(diff_half1_ql,2)*area1 + pow(diff_half2_ql,2)*area2;

    }}

    cout << "Phi separation eval : " << sqrt(sum_errorqs) << "  " << sqrt(sum_errorql) << endl;

    // ========================================================================

    // The average value used to shift pressure
    sumqs /= (double)(mi.MPIglobalCellSize[1] * mi.MPIglobalCellSize[0]);

    double errorqs = 0.0;
    double errorql = 0.0;
    double errorq  = 0.0;
    double errordiff = 0.0;

    // Shift pressure
    for (int j=0; j<mi.MPIglobalCellSize[1]; j++){
    //for (int j=0; j<startj; j++){
    for (int i=0; i<mi.MPIglobalCellSize[0]; i++){

        indice globalcell {i,j};
        int nelem = FlatIndic(mi, globalcell);
        double q, tildeqf;

        fprintf(fqs, "%.16f ", qs_vec.at(nelem));
        fprintf(fqf, "%.16f ", ql_vec.at(nelem));

		  if (j < halfsize){
        qs_exact.at(nelem) = qs_exact.at(nelem) - sum_exactqs_half1 +sum_vecqs_half1;
        } else {
        qs_exact.at(nelem) = qs_exact.at(nelem) - sum_exactqs_half2 +sum_vecqs_half2;
		  }

        fprintf(exactqs, "%.16f ", qs_exact.at(nelem));

        if (mask.at(nelem)){
            ql_exact.at(nelem) = ql_exact.at(nelem) -sum_exactqs_half2 + sum_vecqs_half2;
        }

        fprintf(exactql, "%.16f ", ql_exact.at(nelem));

        // Get Porosity
        vertex local {0.0,0.0};

        basis_->GetCorners(mi, globalcell);

        vertex global = GaussMapPointsFace(local, basis_->corners());

        // Pointwise porosity evaluated at the middle point
        double phif = AssignPorosity(global, pp);
      
        // Averaged porosity of the cell
        double avephi = averagePhi(pp, i, j, mi);

        double tmpq_bar = ql_exact.at(nelem)*avephi + (1.0-avephi)*qs_exact.at(nelem);
        // Accumulate mid-point rules
        double area = mi.cellArea.at(nelem);

        errorqs += pow(qs_vec.at(nelem) - qs_exact.at(nelem), 2)*area;
        errorql += pow(ql_vec.at(nelem) - ql_exact.at(nelem), 2)*area;
        errorq += pow(q_vec.at(nelem) - tmpq_bar, 2)*area;

        if (j<startj){
        double diff = abs(qs_vec.at(nelem) - ql_vec.at(nelem));
        double diffexact = abs(qs_exact.at(nelem) - ql_exact.at(nelem));
      
        errordiff += pow(diff - diffexact,2)* area;}

        if (j>halfsize - cutoff){
        fprintf(errorqs_all, "%.16f ", qs_vec.at(nelem) - qs_exact.at(nelem));
        fprintf(errorql_all, "%.16f ", ql_vec.at(nelem) - ql_exact.at(nelem));
        } else {
        fprintf(errorqs_all, "%.16f ", 0.0);
        fprintf(errorql_all, "%.16f ", 0.0);
        }
    }} 

    cout << "midpoint Error of qs : " << sqrt(errorqs) << endl;
    cout << "midpoint Error of ql : " << sqrt(errorql) << endl;
    cout << "midpoint Error of q : " << sqrt(errorq) << endl;
    cout << "midpoint Error of difference : " << sqrt(errordiff) << endl;


    fclose(exactql);
    fclose(exactqs);
    fclose(fqf);
    fclose(fqs);

    fclose(errorqs_all);
    fclose(errorql_all);

    return 1;
}

/*
int Driver::projectVel(int mark, PhysProperty * pp){

    FILE * dvx = fopen(GetFilename("darcyprox", mark),"w");
    FILE * dvy = fopen(GetFilename("darcyproy", mark),"w");

    for (int j=0; j<mi.MPIglobalCellSize[1]; j++){
    for (int i=0; i<mi.MPIglobalCellSize[0]; i++){

        indice gcell {i,j};
        indice gcellout;
        vertexSet corners = extractCorners(mi, gcell);

    }}

    return 1;
}
*/

// Two sided
int Driver::computeEffVel_case(const vector<vertex>& gaussp,
                               const vertexSet& edgep,
                               const indice& gcellin, const indice& gcellout,
                               const Tensor<weights>& allwgts, double ** lphi,
                               vector<vertex>& vel){
 
    // Extract velocity on these given gauss points
    vector<vertex> vel_relative = 
    ExtractVelocity(&sresult_->vel_darcy, &sresult_->g_darcy,
                    refArrayDarcyEssen_,mi,
                    gaussp, gcellin,*hdiv_,*basis_,{1});
    
    vector<vertex> vel_stokes = 
    ExtractVelocity(&sresult_->vel_stokes, &sresult_->g_stokes,
                    refArrayStokesEssen_,mi,
                    gaussp, gcellin,*br_,*basis_,{1});

    for (int g=0; g<gaussp.size(); g++){
        double phi_mean = 0.0;

        double phiin = advection.eval(gaussp.at(g), ml, location(mi,gcellin), 
                       allwgts({gcellin[0], gcellin[1]}), gcellin, lphi);

        double phiout = advection.eval(gaussp.at(g), ml, location(mi,gcellout), 
                        allwgts({gcellout[0], gcellout[1]}), gcellout, lphi);

        if (phiin == 0.0 && phiout == 0.0){
            phi_mean = 0.0;
        } else {
            phi_mean  = harmonic_mean(phiin ,phiout);
        }

        vel.at(g) = phi_mean*(vel_relative.at(g) + vel_stokes.at(g));

    }

    return 1;
}

// One sided
int Driver::computeEffVel_case(const vector<vertex>& gaussp,
                               const vertexSet& edgep,
                               const indice& gcell,
                               const Tensor<weights>& allwgts, double ** lphi,
                               vector<vertex>& vel){
 
    // Extract velocity on these given gauss points
    vector<vertex> vel_relative = 
    ExtractVelocity(&sresult_->vel_darcy, &sresult_->g_darcy,
                    refArrayDarcyEssen_,mi,
                    gaussp, gcell,*hdiv_,*basis_,{1});
    
    vector<vertex> vel_stokes = 
    ExtractVelocity(&sresult_->vel_stokes, &sresult_->g_stokes,
                    refArrayStokesEssen_,mi,
                    gaussp, gcell,*br_,*basis_,{1});

    for (int g=0; g<gaussp.size(); g++){

        double phi = advection.eval(gaussp.at(g), ml, location(mi,gcell), 
                     allwgts({gcell[0], gcell[1]}), gcell, lphi);

        vel.at(g) = phi*(vel_relative.at(g) + vel_stokes.at(g));

    }

    return 1;
}

int Driver::updateEdgeFlux_case(Tensor<double>& vertedge, Tensor<double>& horiedge,
                                const Tensor<weights>& allwgts, double ** lphi){

    Tensor_zero(vertedge);
    Tensor_zero(horiedge);

    const valarray<double>& gwe = GaussWeightsEdge;
    const valarray<double>& gpe = GaussPointsEdge;

    std::vector<vertex> gaussp;
    gaussp.resize(gpe.size());

    // Effective velocity should be used for transport of concentration
    vector<vertex> effvel; effvel.resize(gaussp.size());

    vector<double> uin;  uin.resize(gaussp.size());
    vector<double> uout; uout.resize(gaussp.size());
 
    for (int j=0; j<mi.MPIglobalCellSize[1]; j++){
    for (int i=0; i<mi.MPIglobalCellSize[0]; i++){

        double flux = 0.0;

        indice gcell {i,j};
        indice cellout;

        effvel.clear(); effvel.resize(gaussp.size()); 

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
           
           computeEffVel_case(gaussp, hori, gcell, allwgts, lphi, effvel);

           flux = edgefluxintegral(hori, CDbottom ,effvel);

        } else {

           cellout = gcell + mi.faceNormal[0];
 
           computeEffVel_case(gaussp, hori, gcell, cellout, allwgts, lphi, effvel);

           flux = edgefluxintegral(mi, gcell, cellout, hori, allwgts, effvel, ml, advection, lphi);

        }        

        horiedge({i,j}) = flux;

        // 1D problem
        vertedge({i,j}) = 0.0;   

    }}

    for (int i=0; i<mi.MPIglobalCellSize[0]; i++){

        // Regarded as outside cell
        indice gcell {i, mi.MPIglobalCellSize[1]-1};
        vertexSet corners = extractCorners(mi, gcell);
        vertexSet hori    = {corners.at(3), corners.at(2)};

        effvel.clear(); effvel.resize(gaussp.size()); 

        // Extract velocity on this edge 
        for (int g=0; g<gpe.size(); g++){
            gaussp.at(g) = GaussMapPointsEdge({gpe[g]},hori);
        }   

        computeEffVel_case(gaussp, hori, gcell, allwgts, lphi, effvel);

        double flux = edgefluxintegral(mi, gcell, hori, allwgts, effvel, ml, advection, lphi);

        if (flux > 0){
            flux = 0.0;
        }

        horiedge({i, mi.MPIglobalCellSize[1]}) = flux; 
    }

    return 1;
}

double melting(const MeshInfo& mi, 
               const vertex& mapped){

    double rate = 0.0;

    if (mapped[1] > -0.2){

        rate = 0.00005*pow(mapped[1]+0.2,1);

    } else {
        rate = 0.0;
    }

    return rate;
}

int Driver::updateCellFlux_case(Tensor<double>& faceflux, 
                                const Tensor<weights>& allwgts, double ** lphi){

    Tensor_zero(faceflux);

    const valarray<double>& gwf = GaussWeightsFace;
    const vector<vertex>& gpf = GaussPointsFace;

    std::vector<vertex> gaussp;
    gaussp.resize(gwf.size());

    for (int j=0; j<mi.MPIglobalCellSize[1]; j++){
    for (int i=0; i<mi.MPIglobalCellSize[0]; i++){

        indice gcell {i,j};
 
        vertexSet corners = extractCorners(mi, gcell); 

        double work = 0.0;

        for (unsigned int g=0; g<gwf.size(); g++){

            double jac = abs(GaussJacobian(gpf[g],corners));
            double gw = gwf[g];
   
            vertex mapped = GaussMapPointsFace(gpf[g],corners);

            work += melting(mi, mapped)* jac * gw;

        }

        double area = mi.cellArea.at(FlatIndic(mi,gcell));

        faceflux({i,j}) = work / area;
    }}

    return 1;
}

int Driver::getflux_case(const Tensor<weights>& allwgts,
                         double ** lphi, double ** lfphi){

    Tensor<double> horiedgeflux = Tensor<double>(2);
    horiedgeflux.setSize({mi.MPIlocalCellSize[0], mi.MPIlocalCellSize[1]+1});

    Tensor<double> vertedgeflux = Tensor<double>(2);
    vertedgeflux.setSize({mi.MPIlocalCellSize[0]+1, mi.MPIlocalCellSize[1]});

    updateEdgeFlux_case(vertedgeflux, horiedgeflux, allwgts, lphi);

    Tensor<double> faceflux = Tensor<double>(2);
    faceflux.setSize({mi.MPIlocalCellSize[0], mi.MPIlocalCellSize[1]});

    updateCellFlux_case(faceflux, allwgts, lphi);

    for (int j=0; j<mi.MPIglobalCellSize[1]; j++){
    for (int i=0; i<mi.MPIglobalCellSize[0]; i++){

        lfphi[j][i] = getcellflux(mi, {i,j}, vertedgeflux, horiedgeflux) - faceflux({i,j});

    }}

    return 1;
}

// ===== Assembly ====
int Driver::CellAvePorosity_case(const indice& gcell, 
                                 const Tensor<weights>& allwgts,
                                 double ** lphi){

    const valarray<double>& gwf = GaussWeightsFace;
    const vector<vertex>&   gpf = GaussPointsFace;

    double phi_f_hat = 0.0;
    double area = 0.0;

    //std::string cellLoc = location(mi, gcell);

    for (unsigned int g=0; g<gwf.size(); g++){
        vertex mapped = GaussMapPointsFace(gpf[g],basis_->corners());

        // Get HD and CD from reconstruction at this gaussian point
        double phif = 0.0;
        phif = abs(advection.eval(mapped, ml, location(mi,gcell), allwgts({gcell[0], gcell[1]}), gcell, lphi));
if (phif<1e-16){phif = 0.0;}

//cout << phif << "  " ;
        // Test =================================================

        phif = AssignPorosity(mapped, myPhase->pp); 

        // ======================================================

        double jac = abs(GaussJacobian(gpf[g],basis_->corners()));
        double gw = gwf[g];
        phi_f_hat += gw * jac * phif;
        area += gw * jac; 
    }

    phi_f_hat /= area;

    myPhase->pp->phi_f_hat = phi_f_hat;
//cout << phi_f_hat << endl;
    return 1;
}

int Driver::AssignLocMatStokes_case(const indice& gcell,
                                    const Tensor<weights>& allwgts,
                                    double ** lphi,
                                    LocMat * loc){

    // copy gaussian quadrature points
    const valarray<double>& gwe = GaussWeightsEdge;
    const valarray<double>& gpe = GaussPointsEdge;
    const valarray<double>& gwf = GaussWeightsFace;
    const vector<vertex>&   gpf = GaussPointsFace;

    // Cell average fluid porosity
    double phi_f_hat = myPhase->pp->phi_f_hat;
    double phi_f = 0.0;
    double phi_s = 0.0;

    // Clear previous calculation
    loc->A.resize(12*12,0.0);
    loc->B.resize(12,0.0);
    loc->f.resize(12,0.0);

    std::fill(loc->A.begin(), loc->A.end(), 0.0);
    std::fill(loc->B.begin(), loc->B.end(), 0.0);
    std::fill(loc->f.begin(), loc->f.end(), 0.0);
    loc->C = 0.0;

    //phi_f_hat = (phi_f_hat == 0.0 ? 1.0 : phi_f_hat);

    for (unsigned int g=0; g<gwf.size(); g++){
        // Calculate mapped gauss points and jacobian
        vertex mapped = GaussMapPointsFace(gpf[g],basis_->corners());
        double jac = abs(GaussJacobian(gpf[g],basis_->corners()));
        double gw = gwf[g];

        // Reconstruction of point wise value of HD and CD
        double phi_f = abs(advection.eval(mapped, ml, location(mi,gcell), allwgts({gcell[0], gcell[1]}), gcell, lphi));
if (phi_f < 1e-16) {phi_f = 0.0;}
//cout << phi_f << "  ";
        // Test ==================================================================

        phi_f = AssignPorosity(mapped, myPhase->pp);
//cout << phi_f << endl;
        // =======================================================================

        phi_s = AssignPorosity(phi_f);       // Solid porosity

        std::array<std::array<double,4>, 12> brwork = 
                           br_->ComputeGradBRmixed(*basis_, mapped);

        std::array<vertex, 12> brval = br_->ComputeBRmixed(*basis_, mapped);

        vertex stokesforce = stokesForce(mapped,myPhase->pp); 

        for (unsigned int j=0; j<12; j++){
                double div1 = brwork[j][0] + brwork[j][3];
            for (unsigned int i=0; i<12; i++){
                double A1 = brwork[j][0];
                double B1 = 0.5*(brwork[j][1] + brwork[j][2]);
                double C1 = brwork[j][3];

                double A2 = brwork[i][0];
                double B2 = 0.5*(brwork[i][1] + brwork[i][2]);
                double C2 = brwork[i][3];
 
                double div2 = brwork[i][0] + brwork[i][3];

                // Symmetrical formulation of A matrix
                loc->A[i+j*12] += 2*phi_s*gw*jac*
                                 (A1*A2+B1*B2*2+C1*C2 - (1.0/3.0)*div1*div2);
            }
            // With dimension version
            loc->B[j] += gw*jac*div1 * br_->Pressure();

            // Non dimensionalized version
            // Attention, porosity has been multiplied to right hand side force term
            loc->f[j] += gw*jac* phi_s*(stokesforce[0]*brval[j][0] + 
                                        stokesforce[1]*brval[j][1]);

        }

        // Non dimensionalized version
//        loc->C += gw*jac*phi_f_hat/phi_s*
//                  br_->Pressure()*br_->Pressure();
        loc->C += gw*jac*phi_f/phi_s*
                  br_->Pressure()*br_->Pressure();
    }

    return 1;
}

inline bool outside(const MeshInfo& mi, const indice& cell){

    if (cell[0] < 0 || cell[1] < 0 || cell[1] > mi.MPIglobalCellSize[1]-1 || cell[0] > mi.MPIglobalCellSize[0]-1){
        return true;
    } else {
        return false;
    }
}

int Driver::AssignLocMatDarcy_case(const indice& gcell,
                                   const Tensor<weights>& allwgts,
                                   double ** lphi,
                                   LocMat * loc){
                     
    // copy gaussian quadrature points
    const valarray<double>& gwe = GaussWeightsEdge;
    const valarray<double>& gpe = GaussPointsEdge;
    const valarray<double>& gwf = GaussWeightsFace;
    const vector<vertex>&   gpf = GaussPointsFace;

    // Cell average fluid porosity
    double phi_f_hat = myPhase->pp->phi_f_hat;
    double phi_f = 0.0;
    double phi_s = 0.0;
    double theta = myPhase->pp->theta;

    loc->A.resize(8*8,0.0);
    loc->B.resize(8,0.0);
    loc->f.resize(8,0.0);

    std::fill(loc->A.begin(), loc->A.end(), 0.0);
    std::fill(loc->B.begin(), loc->B.end(), 0.0);
    std::fill(loc->f.begin(), loc->f.end(), 0.0);
    loc->C = 0.0;

    for (unsigned int g=0; g<gwf.size(); g++){
        // Calculate mapped gauss points and jacobian
        vertex mapped = GaussMapPointsFace(gpf[g],basis_->corners());
        double jac = abs(GaussJacobian(gpf[g],basis_->corners()));
        double gw = gwf[g];

        // Reconstruction of point wise value of HD and CD
        double phi_f = abs(advection.eval(mapped, ml, location(mi, gcell), allwgts({gcell[0], gcell[1]}), gcell, lphi));
if (phi_f < 1e-16) {phi_f = 0.0;}
//cout << phi_f << "  " ;
        // Test ==================================================================

        phi_f = AssignPorosity(mapped, myPhase->pp);
//cout << phi_f << endl;
		  // =======================================================================

        phi_s = AssignPorosity(phi_f);

        std::array<vertex, 8> hdivwork = hdiv_->ComputeHdivmixed(*basis_,mapped);

        vertex darcyforce = darcyForce(mapped,myPhase->pp);

        for (unsigned int j=0; j<8; j++){
            for (unsigned int i=0; i<8; i++){

                // Non dimensionalized version
                loc->A[i+j*8] += gw*jac* 
                                 (hdivwork[j][0]*hdivwork[i][0] + 
                                  hdivwork[j][1]*hdivwork[i][1]);
            }
            // darctforce is set to be zero here
            loc->f[j] += gw*jac*(darcyforce[0]*hdivwork[j][0] + 
                                 darcyforce[1]*hdivwork[j][1]);
        }

        // Compaction matrix
        //loc->C += gw*jac*1.0/phi_s*
        //          hdiv_->Pressure()*hdiv_->Pressure();

        double scaletmp = 1.0;
        if (phi_f_hat > 1e-15){
            scaletmp = phi_f/phi_f_hat;
        }

        loc->C += gw*jac*scaletmp/phi_s*
                  hdiv_->Pressure()*hdiv_->Pressure();

    }

/*
    if (gcell[1] == 34 || gcell[1] == 35){
phi_f_hat = 2.0/30.0;
    }
cout << phi_f_hat << endl;
*/
    // Define B matrix for the Darcy part
    // Compute with divergence theorem
    phi_f_hat = (phi_f_hat == 0.0 ? 1.0 : phi_f_hat);

    vertexSet corners = basis_->corners();

    for (int e =0; e<4; e++){

        vertexSet corner = {corners.at((e+3)%4),
                            corners.at(e)};
        double len = length(corner);
        for (int g=0; g<gpe.size(); g++){
            vertex mapped = GaussMapPointsEdge({gpe[g]},corner);
            std::array<vertex, 8>  hdivwork = hdiv_->ComputeHdivmixed(*basis_,mapped);
            // Zeroth order constant pressure basis is always 1
            vertex nu = basis_->unitnormal(e);

            indice cellout = gcell + mi.faceNormal[((e-1)+4)%4];

            double phi_f_e = 0.0;
/*
            if (outside(mi, cellout)){
                phi_f_e = abs(advection.eval(mapped, ml, location(mi, gcell), allwgts({gcell[0], gcell[1]}), gcell, lphi));
if (phi_f_e < 1e-16) {phi_f_e = 0.0;}

            } else {
                double in = abs(advection.eval(mapped, ml, location(mi, gcell), allwgts({gcell[0], gcell[1]}), gcell, lphi)); 
                double out = abs(advection.eval(mapped, ml, location(mi, cellout), 
                             allwgts({cellout[0], cellout[1]}), cellout, lphi)); 
					 if (in < 1e-16){in = 0.0;}
					 if (out< 1e-16){out = 0.0;}
//cout<< endl << in << "  " << out << endl;
        if (in < 1e-16 && out <1e-16){
            phi_f_e = 0.0;
        } else {
            phi_f_e  = harmonic_mean(in ,out);
        }

            }
*/

            // Reconstruction of point wise value of HD and CD
//cout << gcell[0] << "  " << gcell[1] << "  " << phi_f_e << "  ";
            // Testing =================================================

            phi_f_e = AssignPorosity(mapped, myPhase->pp);

/*
				if ((gcell[1]==35 && e==1) ||(gcell[1]==34 && e==3)){
            phi_f_e = 2.0/(1.0/0.1 + 1.0/0.05);
            } 
*/
//cout << gcell[0] << "  " << gcell[1] << "  " << e << "  " <<  phi_f_e << "  " << endl;
            // =========================================================

            for (int j=0; j<8; j++){
                // With dimension version
                loc->B[j] += len/2.0*gwe[g]*
                             pow(phi_f_hat,-0.5) * pow(phi_f_e, 1+theta) *
                             (hdivwork[j][0] * nu[0]+
                              hdivwork[j][1] * nu[1]);
            } 
        }
    }

    return 0;
}

int Driver::AssignLocMatCouple_case(const indice& gcell,
                                    const Tensor<weights>& allwgts,
                                    double ** lphi,
                                    double& k){

    const valarray<double>& gwf = GaussWeightsFace;
    const vector<vertex>&   gpf = GaussPointsFace;

    k = 0.0;

    double phi_f_hat = myPhase->pp->phi_f_hat;
    double phi_f = 0.0;
    double phi_s = 0.0;

    for (unsigned int g=0; g<gwf.size(); g++){
        // Calculate mapped gauss points and jacobian
        vertex mapped = GaussMapPointsFace(gpf[g],basis_->corners());
        double jac = abs(GaussJacobian(gpf[g],basis_->corners()));
        double gw = gwf[g];

        // Reconstruction of point wise value of HD and CD
        double phi_f = abs(advection.eval(mapped, ml, location(mi, gcell), allwgts({gcell[0], gcell[1]}), gcell, lphi));
if (phi_f < 1e-16) {phi_f = 0.0;}
//cout << phi_f << "  ";
        // Test ============================================================

        phi_f = AssignPorosity(mapped, myPhase->pp);
//cout << phi_f << endl;
        // =================================================================

        phi_s = AssignPorosity(phi_f);

        //k -= gw*jac*pow(phi_f_hat,0.5)/phi_s * br_->Pressure() * 
        //                                       hdiv_->Pressure();

        double scaletmp = 0.0;
        if (phi_f_hat > 1e-16){
            scaletmp = phi_f/sqrt(phi_f_hat);
        }

        k -= gw*jac*scaletmp/phi_s * br_->Pressure() * 
                                     hdiv_->Pressure();

    }

    return 1;
}

int Driver::PrepareTransport_case(double (*func)(const valarray<double>& point,
                                                 const vector<double>& param) ){

    PetscCall(DMCreateGlobalVector(dmu, &globalCD));

    SimpleInitialValue(dmMesh, dmu, &globalmesh, &globalCD, {H_,0.0}, func);

    // Initialization of multi level weno and corresponding usage
    ml = multilevel(); 

    ml.addLevel("(1,3)", {1,3}, mi);
    ml.addLevel("(1,2)", {1,2}, mi);

    // Area scale
    h0 = sqrt((L_*H_)/
         (double)(mi.MPIglobalCellSize[0]*mi.MPIglobalCellSize[1]));

    advection = mluse();

    unordered_map<std::string, vector<indice>> method;
    method.insert(std::make_pair<std::string, vector<indice>>("(1,3)", { {0,-1} }));
    method.insert(std::make_pair<std::string, vector<indice>>("(1,2)", { {0,-1} , {0,0} }));

    advection.setmethod("all", method);
    advection.setbias("all");

    CDbottom = func({0.0,-1*H_}, {H_,0.0});

    return 1;
}

int Driver::RK_case(double dt, double Tmax, int maxIter, double tolUzawa){

    int mark = 1 + start;

    int Nt = (int)(Tmax/dt);

    for (int t=0; t<Nt; t++) {

        Vec localphi;
        PetscCall(DMGetLocalVector(dmu, &localphi));

        Vec fluxphi;
        PetscCall(VecDuplicate(globalCD, &fluxphi));
 
        double ** lphi;
        double ** lfphi;

        PetscCall(DMGlobalToLocalBegin(dmu, globalCD, INSERT_VALUES, localphi));
        PetscCall(DMGlobalToLocalEnd(dmu, globalCD, INSERT_VALUES, localphi));

        PetscCall(DMDAVecGetArray(dmu, localphi, &lphi););
        PetscCall(DMDAVecGetArray(dmu, fluxphi, &lfphi));

        ml.updatesigma(lphi);
        Tensor<weights> allwgts;
        advection.computeWgts(ml, mi, h0, allwgts, location);

        cout << "Darcy-Stokes system solved at : " << t*dt << endl;
        SolveFlow_case(maxIter, tolUzawa, allwgts, lphi);
        CreateScatterVec();

        getflux_case(allwgts, lphi, lfphi);

        DMDAVecRestoreArray(dmu, fluxphi, &lfphi);
        DMDAVecRestoreArray(dmu, localphi, &lphi);
        DMRestoreLocalVector(dmu, &localphi);

        //printCellAve(mark, &globalCD, mi, "porosity");
        printExactPorosity(mark, mi, "porosity", myPhase->pp);
        printVelEdgeGauss_case(mark, myPhase->pp);

        //printSimplePressure_case(mark, myPhase->pp);
        printShiftedPressure_case(mark, myPhase->pp);
        mark ++;
        cout << errorNorm(mark, myPhase->pp, 2) << endl;

//VecView(fluxphi, PETSC_VIEWER_STDOUT_WORLD);
        VecAXPY(globalCD, -1*dt, fluxphi);
    }

    return 1;
}

int Driver::ParallelMatrixAssemble_case(const Tensor<weights>& allwgts,
                                        double ** lphi){

    PetscMPIInt size, rank;

    MPI_Comm_size(PETSC_COMM_WORLD, &size);
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    PetscFunctionBeginUser;

    // Get gauss points first
    const valarray<double>& gwe = GaussWeightsEdge;
    const valarray<double>& gpe = GaussPointsEdge;
    const valarray<double>& gwf = GaussWeightsFace;
    const vector<vertex>&   gpf = GaussPointsFace;

    // Calculate dofs 
    int totalElem = mi.MPIglobalCellSize[0] * mi.MPIglobalCellSize[1];

    int reducedDOFStokes = br_->getDOF() - bndryDOFStokes_;

    int reducedDOFDarcy = hdiv_->getDOF() - bndryDOFDarcy_;

    PrepareReducedSys(reducedStokes_, reducedDOFStokes, bndryDOFStokes_, 
                      totalElem, 30, 30, 4, 4);
    PrepareReducedSys(reducedDarcy_, reducedDOFDarcy, bndryDOFDarcy_, 
                      totalElem, 14, 14, 2, 2);

    PetscCall(MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, 
                           totalElem, totalElem, 
                           1, NULL, 0, NULL, &K));
    PetscCall(MatSetUp(K));

    // ===================================================================

    LocMat * locmatS = new LocMat;
    LocMat * locmatD = new LocMat;

    double k = 0.0;

    // ! Loop local portion of physical domain
    int istart = mi.MPIlocalCellStart[0];
    int jstart = mi.MPIlocalCellStart[1];

    for (int j=jstart; j<jstart + mi.MPIlocalCellSize[1]; j++){
    for (int i=istart; i<istart + mi.MPIlocalCellSize[0]; i++){

        // ! Get global element index
        indice global {i,j};

        int nElem = FlatIndic(mi,global);

        // ! Extract corners of this element
        basis_->GetCorners(mi, global);

        // ! Compute cell averaged porosity
        CellAvePorosity_case(global, allwgts, lphi);

        // ! Compute local values associated to each dofs
        AssignLocMatStokes_case(global, allwgts, lphi, locmatS);
        AssignLocMatDarcy_case(global, allwgts, lphi,locmatD);
        AssignLocMatCouple_case(global, allwgts, lphi, k);

        // ! Load corresponding shape functions
        shape stokesFuncSp(basis_, br_);
        shape darcyFuncSp(basis_, hdiv_);

        // ! Assign local values to global matrix
        if (elemOnBndry(mi, global)){
            // ! Dealing wiht boundary dofs
            AssignLocRedSys(reducedStokes_, locmatS, refArrayStokesEssen_, 
                            mi, bndryStokesEssen_, global, stokesFuncSp, parameter); 
            AssignLocRedSys(reducedDarcy_, locmatD, refArrayDarcyEssen_,
                            mi, bndryDarcyEssen_, global, darcyFuncSp, parameter);
        } else {
            AssignLocRedSys(reducedStokes_, locmatS, refArrayStokesEssen_, mi, global, *br_);
            AssignLocRedSys(reducedDarcy_, locmatD, refArrayDarcyEssen_, mi, global, *hdiv_);
        }

        // Assign coupling K matrix and two C matrices
        // const pressure space not affected by boundary dofs
        PetscCall(MatSetValue(K,nElem,nElem,k,ADD_VALUES));
        PetscCall(MatSetValue(reducedStokes_->C, nElem, nElem, locmatS->C,ADD_VALUES));
        PetscCall(MatSetValue(reducedDarcy_->C, nElem, nElem, locmatD->C, ADD_VALUES));

    }}

    AssembleReducedSys(reducedStokes_);
    AssembleReducedSys(reducedDarcy_);

    PetscCall(MatAssemblyBegin(K,MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(K,MAT_FINAL_ASSEMBLY));

    return 1;
}

int Driver::SolveFlow_case(int maxIter, double tolUzawa, 
                           const Tensor<weights>& allwgts, double ** lphi){

    ParallelMatrixAssemble_case(allwgts, lphi);

    int nelem = mi.MPIglobalCellSize[0] * mi.MPIglobalCellSize[1];

    CreateLinearSys(reducedStokes_, nelem);
    CreateLinearSys(reducedDarcy_, nelem);

    CreateCoupledSystem(reducedStokes_, reducedDarcy_, Result_, &K);

    CoupledUzawa(Result_, tolUzawa, maxIter);

    // New exact solver
    //SchurSolver(Result_);

    return 1;
}
