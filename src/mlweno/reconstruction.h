#ifndef RECONSTMLWENO_H_
#define RECONSTMLWENO_H_

#include "stencilpolynomial.h"
#include <unordered_map>

class reconstruction{

    public:
        reconstruction() {};
        ~reconstruction() {};

        int prepare(const vector<int>& insize, const MeshInfo& mi);

        double eval(const vector<int>& index, 
                    const Tensor<double>& sol,
                    const vertex& point) const;

        double eval(const vector<int>& stencilindex,
                    const vector<int>& baseindex,
                    const vertex& point) const;

        double sigma(const vector<int>& index,
                     const Tensor<double>& sol) const;

        int dsigma(const vector<int>& index,
                   const Tensor<double>& sol,
                   vector<double>& der) const;

        int getSize() const{return stencilPoly.getSize();};
        int getSize(const int& dim) const{return stencilPoly.getSize(dim);};

        int getorder() const{return stencilPoly(0).getorder();};

        int getStencilSize(const int& dim) const{return size.at(dim);};

        int printcoef(){for (int i=0;i<stencilPoly.getSize(); i++){stencilPoly(i).printCoef();} return 1;}

        int printCoef(int s){return stencilPoly(s).printCoef();}

        int printSigmaTensor(int s);

    private:

        // Reconstruction stencil size
        vector<int> size {-1,-1};

        Tensor<stencilpolynomial> stencilPoly;

        int getcenter(const vector<vector<vertex>>& cornerSet,
                      vector<vertex>& refcell,
                      vertex& center, double& h, double& area);

};

class multilevel{

    public:
        multilevel() {};
        ~multilevel() {};

        int addLevel(const std::string& name, 
                     const vector<int>& stencilSize,
                     const MeshInfo& mi);

        double eval(const std::string& name, 
                    const vector<int>& index, 
                    const Tensor<double>& sol,
                    const vertex& point) const;

        double eval(const std::string& name, 
                    const vector<int>& stencilindex,
                    const vector<int>& baseindex,
                    const vertex& point) const;

        double sigma(const std::string& name, 
                     const vector<int>& index,
                     const Tensor<double>& sol) const;

        int dsigma(const std::string& name, 
                   const vector<int>& index,
                   const Tensor<double>& sol,
                   vector<double>& der) const;

        int getorder(const std::string& name) const
        {return mlrecons.at(name).getorder();}

        int getSize(const std::string& name)const 
        {return mlrecons.at(name).getSize();};

        int getSize(const std::string& name, const int& index)const 
        {return mlrecons.at(name).getSize(index);};

        int getStencilSize(const std::string& name, const int& index)const 
        {return mlrecons.at(name).getStencilSize(index);};

        int printSigmaTensor(const std::string& name, const int& index)
        {return mlrecons.at(name).printSigmaTensor(index);};

        int printCoef(const std::string& name, const int& index)
        {return mlrecons.at(name).printCoef(index);};

        /**!
         * Extract cell averaged solution to a tensor object
         */
        int getsol(Tensor<double>& stencilsol, double ** localsol,
                   const indice& stencilindex, const std::string& name) const;

        /**!
         * Update smoothness indicator of all levels with cell-averaged value
         */
        int updatesigma(double ** localsol);

        double getsigma(const std::string& name, const vector<int>& index) const
        {return alllevelsigma.at(name)(index);};

        double getscaledsigma(const std::string& name, const vector<int>& index)const
        {return alllevelscaled.at(name)(index);}

        int updatedersigma(double ** localsol, const MeshInfo& mi);

        unordered_map<int,double> getdersigma(const std::string& name, const vector<int>& index) const
        {return allleveldersigma.at(name)(index);}

        unordered_map<int,double> getderscaledsigma(const std::string& name, const vector<int>& index) const
        {return alllevelderscaled.at(name)(index);}

        // Set holds all the reconstruction levels
        set<std::string> reconlevelSet;

        int geteta(const int& rl) const;

        /**!
         * Calculate everything at once including sigma, derivative of sigma,
         * scaled sigma, and derivative of scaled sigma
         */
        int updateall(double ** localsol, const double& h0, const int& s, const double& ep, const MeshInfo& mi) ;

        // =======================================================================================================
        /**!
         * Print out all stencil polynomial coefficients
         */
        int printcoef(const std::string& name);



        int printsigma(const std::string& name);

        int printscaledsigma(const std::string& name);

        int printdsigma(const std::string& name);

        int printdscaledsigma(const std::string& name);

    private:

        unordered_map<std::string, reconstruction> mlrecons;
        unordered_map<std::string, Tensor<double>> alllevelsigma;
        unordered_map<std::string, Tensor<unordered_map<int,double>>> allleveldersigma;

        // Compute identifier utilizing smoothness indicator and derivative of smoothness indicator
        unordered_map<std::string, Tensor<double>> alllevelscaled;
        unordered_map<std::string, Tensor<unordered_map<int,double>>>  alllevelderscaled;
};

#endif
