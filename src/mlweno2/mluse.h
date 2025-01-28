#ifndef MLUSE_H_
#define MLUSE_H_
/**!
 * For each usage of multi-level weno reconstruction
 * one mluse object is needed.
 * For example, for advection (3,2,1) is used and considered a mluse object.
 * Sigam, jacobian and non-linear wgts are all stored in this object.
 */
#include "reconstruction.h"

class mluse {

    public:
        mluse() { ep = 1e-4; s = 1;};
        ~mluse() {};

        /**!
         * Set reconstruction method.
         */
        int setmethod(const std::string& pos,
                      const unordered_map<std::string, vector<indice>>& method);

        /**!
         * Compute non linear weights
         * with mlweno weighting scheme 
         * should be called after setmethod being called
         */

        int setbias(const std::string& pos);


        //int getsol(Tensor<double>& stencilsol, double** localsol, 
        //           const indice& stencilindex) const;

        /**!
         * Update smoothness indicator for all recorded levels.
         */
        //int updatesigma(const multilevel& ml, double ** localsol);

        int printsigma(const std::string& name);

//        int computeWgts(const unordered_map<std::string, vector<indice>>& method,
//                        unordered_map<std::string, vector<double>>& wgts, 
//                        const multilevel& ml, const double& h0,
//                        const indice& index);

        int computeWgts(const std::string& pos, const multilevel& ml,
                        const indice& index,    const double& h0,
                        unordered_map<std::string, vector<double>>& wgts);

        int printWgts(const unordered_map<std::string, vector<double>>& wgts);

        /**!
         * Evaluation of the nonllinear weighted value at the given point.
         */
        //double eval(const vertex& point, const multilevel& ml,
        //            const unordered_map<std::string, vector<indice>>& method,
        //            const unordered_map<std::string, vector<double>>& wgts,
        //            const indice& index, double ** localsol)const;

        double eval(const vertex& point, const multilevel& ml,
                    const std::string& pos, 
                    const unordered_map<std::string, vector<double>>& wgts,
                    const indice& index, double ** localsol) const;

    public:

        set<std::string> posSet;

        unordered_map<std::string, unordered_map<std::string, vector<indice>>> reconstMethod;

        unordered_map<std::string, Tensor<double>> sigma;

        double ep = 1e-4;
        int    s  = 1;

        unordered_map<std::string, unordered_map<std::string, double>> bias;

        int geteta (const int& rl) const;

        bool stencilexist(const multilevel& ml, const indice& index, const std::string& name) const;
};

#endif
