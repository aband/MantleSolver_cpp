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

        double sigma(const vector<int>& index,
                     const Tensor<double>& sol) const;

        int getSize() const{return stencilPoly.getSize();};
        int getSize(const int& dim) const{return stencilPoly.getSize(dim);};

        int getStencilSize(const int& dim) const{return size.at(dim);};

        int printcoef(){for (int i=0;i<stencilPoly.getSize(); i++){stencilPoly(i).printCoef();} return 1;}

        // Derivative

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

        double sigma(const std::string& name, 
                     const vector<int>& index,
                     const Tensor<double>& sol) const;

        int getSize(const std::string& name)const 
        {return mlrecons.at(name).getSize();};

        int getSize(const std::string& name, const int& index)const 
        {return mlrecons.at(name).getSize(index);};

        int getStencilSize(const std::string& name, const int& index)const 
        {return mlrecons.at(name).getStencilSize(index);};

        /**!
         * Extract cell averaged solution to a tensor object
         */
        int getsol(Tensor<double>& stencilsol, double ** localsol,
                   const indice& stencilindex, const std::string& name) const;

        /**!
         * Update smoothness indicator of all levels with cell-averaged value
         */
        int updatesigma(double ** localsol);

        int printsigma(const std::string& name);

        double getsigma(const std::string& name, const vector<int>& index) const
        {return alllevelsigma.at(name)(index);};

        /**!
         * Print out all stencil polynomial coefficients
         */
        int printcoef(const std::string& name);

        // Set holds all the reconstruction levels
        set<std::string> reconlevelSet;

    private:

        unordered_map<std::string, reconstruction> mlrecons;
        unordered_map<std::string, Tensor<double>> alllevelsigma;
};

#endif
