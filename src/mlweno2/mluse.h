#ifndef MLUSE_H_
#define MLUSE_H_
/**!
 * For each usage of multi-level weno reconstruction
 * one mluse object is needed.
 * For example, for advection (3,2,1) is used and considered a mluse object.
 * Sigam, jacobian and non-linear wgts are all stored in this object.
 */
#include "reconstruction.h"

using weights = std::unordered_map<std::string, std::vector<double>>;

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
         * Set linear weights
         * with mlweno weighting scheme 
         * should be called after setmethod being called
         */
        int setbias(const std::string& pos);

        /**!
         * Compute non linear weighting
         */
        int computeWgts(const std::string& pos, const multilevel& ml,
                        const indice& index,    const double& h0,
                        weights& wgts);

        int computeWgts(const std::string& pos, const indice& gcell, 
                        const multilevel& ml, weights& wgts);

        /**!
         * Compute non linear weighting for all stencils at once
         * and store in a tensor object.
         */
        int computeWgts(const multilevel& ml, const MeshInfo& mi, const double& h0, Tensor<weights>& allwgts);

        int computeWgts(const multilevel& ml, const MeshInfo& mi, Tensor<weights>& allwgts);

        /**!
         * Print selected non linear weights.
         */

        int printWgts(const weights& wgts);

        int printWgts(const Tensor<weights>& allwgts, const vector<int>& index) {return printWgts(allwgts(index));};

        /**!
         * Evaluation of the nonllinear weighted value at the given point.
         */
        double eval(const vertex& point, const multilevel& ml,
                    const std::string& pos, 
                    const weights& wgts,
                    const indice& index, double ** localsol) const;

        int sumweights(const multilevel& ml, const mluse& use,
                       double& sumwgts,
                       derivative& sumderwgts,
                       const std::string& pos,
                       const indice& gcell) const;

        int der(const vertex& point, const multilevel& ml,
                const::string& pos,  const weights& wgts, 
                const indice& index, double ** localsol,
                const MeshInfo& mi,
                unordered_map<int, double>& derivative) const;

    public:

        set<std::string> posSet;

        unordered_map<std::string, unordered_map<std::string, vector<indice>>> reconstMethod;

        double ep = 1e-4;
        int    s  = 1;

        unordered_map<std::string, unordered_map<std::string, double>> bias;

        int geteta (const int& rl) const;

        bool stencilexist(const multilevel& ml, const indice& index, const std::string& name) const;
};

#endif
