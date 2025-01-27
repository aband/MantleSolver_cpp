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

        // Set holds all the reconstruction levels
        set<std::string> reconlevelSet;

    private:

        unordered_map<std::string, reconstruction> mlrecons;
};

/**!
 * For each usage of multi-level weno reconstruction
 * one mluse object is needed.
 * For example, for advection (3,2,1) is used and considered a mluse object.
 * Sigam, jacobian and non-linear wgts are all stored in this object.
 */
class mluse {

    public:
        mluse() {};
        ~mluse() {};

        int setmethod(const std::string& pos,
                      const unordered_map<std::string, vector<indice>>& method);

        int getsol(Tensor<double>& stencilsol, double** localsol, 
                   const indice& stencilindex);

        /**!
         * Update smoothness indicator for all recorded levels.
         */
        int updatesigma(const multilevel& ml, double ** localsol);

        int printsigma(const std::string& name);

        /**!
         * Compute non linear weight
         */
        int computeWgts(const unordered_map<std::string, vector<indce>>& method,
                        unordered_map<std::string, vector<double>>& wgts);

    public:

        set<std::string> posSet;

        unordered_map<std::string, unordered_map<std::string, vector<indice>>> reconstMethod;

        unordered_map<std::string, Tensor<double>> sigma;
};

#endif
