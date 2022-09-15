
// Couple of constant functions
double constfunc(valarray<double>& point,const vector<double>& param){
    return 1.0;
}

double constfunc(valarray<double>& point,const vector<double>& param, double c){
    return c;
}

double constfunc(){
    return 1.0;
}

double constfunc(double c){
    return c;
}


// Evaluation of factorial
int factorial(int top, int bottom){
    assert(top>bottom || top==bottom);
    if (top==bottom){ return 1;}
    else{ return top*factorial(top-1,bottom);}
}

int factorial(int top){
    assert(top>0 || top==0);
    if (top==0){ return 1;}
    else{ return top*factorial(top-1);}
}

