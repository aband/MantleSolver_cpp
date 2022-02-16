// ====================================================================================
class WenoBasisCoeff : public WenoMesh{
    public:
        // Constructor
        WenoBasisCoeff(int M, int N, int g, vector <valarray<double> >& lm) :
                       WenoMesh(M, N, g, lm) {};

        void SetStencilInfo(vector< valarray<int> >& stencilindex, 
                            vector<int>& order, vector<int>& range);

        // Get stencil information. Need to be called first!
        void GetStencil();

        // Create weno basis polynomial coefficients based on the given stencil
        void CreateWenoBasisCoeff();

        // Print out calculated weno basis coefficients for testing
        void PrintWenoBasisCoeff();

        // Print out assigned stencil for a weno basis 
        void PrintStencil();

    protected:

        vector< vector< vector<double> > > * wenobasiscoeff;

    private:
        vector<int> order;       // WENO restruction polynomial order, (xpow, ypow)
         /*
         *Index set used in calculation of the reconstruction is orderred in the
         *following. The members of this vector will be added to the starting index.
         */
        vector< valarray<int> > stencilindex;     // Index set used in reconstruction

        /*
         *range of index that include start index and end index of i and j.
         *ibegin,iend,jbegin,jend.
         */
        vector<int> range;

        /*
         *Stencil contains four corners of each cell within the defined
         *stencil. Center point is defined by averaging out all points
         *in stencil.
         */
        vector< vector< vector< valarray<double> > > > * stencil;
        vector< valarray<double> > * stencilcenter;
        vector<double> * stencilh;

};

class SmoothnessIndicator : public WenoMesh{
    public:
        SmoothnessIndicator(int M, int N, int g, vector <valarray<double> >& lm) :
                            WenoMesh(M, N, g, lm) {};
    protected:

};

class WenoReconstruction : public WenoBasisCoeff{
    public:
        WenoReconstruction(int M, int N, int g, vector <valarray<double> >& lm) :
                           WenoBasisCoeff(M, N, g, lm) {};


};

// =============================================================
void WenoBasisCoeff::SetStencilInfo(vector< valarray<int> >& si, 
                    vector<int>& o, vector<int>& r){
    order = o;
    stencilindex = si;
    range = r;
}

void WenoBasisCoeff::GetStencil(){

    stencil = new vector< vector< vector< valarray<double> > > >();
    stencilcenter = new vector< valarray<double> >();
    stencilh = new vector<double>();

    int totali = M+2*ghost;

    for (int j=range[2]; j<range[3]; j++){
    for (int i=range[0]; i<range[1]; i++){

        // Locate reference center point for each stencil
        vector< vector< valarray<double> > > stencilcell;
        valarray<double> scenter {0.0,0.0};

        for (auto & s: stencilindex){
            vector< valarray<double> > cellcorner;
            valarray<double> cellcenter {0.0,0.0};

            for (auto & c: cornerindex){
                valarray<double> corner = lmesh[(j+c[1]+s[1])*totali+(i+c[0]+s[0])];
                cellcenter += corner/(4.0); 
                cellcorner.push_back(corner);
            }

            scenter +=cellcenter/(double)(order[0]*order[1]);

            stencilcell.push_back(cellcorner);
        }
        vector <double> Empty;

        double h = pow(NumIntegralFace(stencilcell.at(0),Empty,{0.0,0.0},1.0,constfunc),0.5);

        stencil->push_back(stencilcell);
        stencilcenter->push_back(scenter);  
        stencilh->push_back(h);
    }}
}

void WenoBasisCoeff::CreateWenoBasisCoeff(){

    wenobasiscoeff = new vector< vector< vector<double> > >();

    for (int i=0; i<stencil->size(); i++){

        vector< vector<double> > stencilcoeff;
        /*
         *Get corresponding reference center point and normalize
         *factor h for each stencil.
         */
        valarray<double> center = stencilcenter->at(i);
        double h = stencilh->at(i);

        /*
         *Preparing parameters for solving linear system
         *using lapack.
         */
        lapack_int n = order[0]*order[1];
        lapack_int nrhs = n; 
        lapack_int lda = n; 
        lapack_int ldb = nrhs;

        double * a = new double [n*n];
        double * b = new double [n*nrhs];
        lapack_int * p = new int [n];

        // Loop through a single stencil
        for (int cell=0; cell<stencil->at(i).size(); cell++){

            for (int ypow=0; ypow<order[1]; ypow++){
            for (int xpow=0; xpow<order[0]; xpow++){
                vector<double> param {(double)xpow,(double)ypow}; 
                a[cell*n + ypow*order[0] + xpow] = NumIntegralFace(stencil->at(i).at(cell),param,center,h,poly);
            }} 

        }

        fill(b,b+n*nrhs,0);

        for (int i =0; i<nrhs; i++){b[i*n+i] = a[n*i];}

        int err = LAPACKE_dgesv(LAPACK_ROW_MAJOR, n, nrhs, a, lda, p, b, ldb);
        if (err){
            printf("ERROR: Weno Basis Coefficient for order %d, %d. Error type %d \n",
                   order[0],order[1],err);
        }
        vector<double> work;
        work.assign(b,b+n*nrhs);
        stencilcoeff.push_back(work);

        wenobasiscoeff->push_back(stencilcoeff);
    }
}

void WenoBasisCoeff::PrintStencil(){

    printf("The total number of stencil is %ld \n",stencil->size());
    printf("The number of cell of each stencil is %ld \n",stencil->at(0).size());

    for (int s=0; s<stencil->size(); s++){
        for (auto & cell: stencil->at(s)){
            for (auto & corner: cell){
                printf("(%.2f,%.2f),  ",corner[0],corner[1]);
            } cout << endl;
        } cout << endl;
    } cout << endl;
}

void WenoBasisCoeff::PrintWenoBasisCoeff(){

    for (int s=0; s<wenobasiscoeff->size(); s++){
         for (auto & cell: wenobasiscoeff->at(s)){
             int n = order[0]*order[1];
             cell.resize(n,n);
             for (int j=0; j<n; j++){
             for (int i=0; i<n; i++){
                 printf("coeff = %.4f ", cell[i*n+j]);
             } cout << endl;}cout<< endl;
         } cout << endl;
    } cout << endl;
}
