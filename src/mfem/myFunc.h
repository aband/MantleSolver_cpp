#ifndef MYFUNC_H_
#define MYFUNC_H_

typedef struct {

    double theta;
    double mu_s;
    double mu_f;
    double rho_f = 2800;
    double rho_s = 3300;
    double gx = 0.0;
    double gy = -10.0;

    double l = 20;
} PhysProperty;

double AssignPorosity(const vertex& point, const double& l){

    if (point[1] < 12000 && abs(point[0]) < point[1] + l){
        return 0.05*pow(1.0-point[1]/120000,2) * (1-abs(point[0])/(l+point[1]));
    } else {
        return 0.0;
    }
}

#endif
