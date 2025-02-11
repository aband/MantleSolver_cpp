#include "advectiveflux.h"

inline double LFflux(const double& uin, const double& uout, 
                     const double& fin, const double& fout,
                     const double& LF){

    return 0.5*(fout + fin - LF*(uout - uin));
}

inline int derLFflux(const derivative& derin,  const derivative& derout,
                     const derivative& fderin, const derivative& fderout,
                     const double& LF,
                     derivative& work){

    // differentiating Lax_Friedrich flux

    work = fderin;
    unordered_map_arithmetic(work, fderout, std::plus<double>());

    unordered_map_arithmetic(work, derout, std::minus<double>(), 
                                   LF, std::multiplies<double>());

    unordered_map_arithmetic(work, derin,  std::plus<double>(), 
                                   LF, std::multiplies<double>());

    unordered_map_arithmetic(work, 0.5, std::multiplies<double>());

    return 1;
}

double edgefluxintegral(const MeshInfo& mi, 
                        const indice& gcellin,
                        const indice& gcellout,
                        const vertexSet& edge,
                        const Tensor<weights>& allwgts,
                        const vector<vertex>& vel,
                        multilevel& ml,
                        mluse& use,
                        double ** lu){

    double work = 0.0;

    //! Extract default gauess points and gauess weights.
    const valarray<double>& gwe = GaussWeightsEdge;
    const valarray<double>& gpe = GaussPointsEdge;

    // Get edge lendth and unit vector normal to the given edge
    double len = length(edge);
    vertex unitNormal = UnitNormal(edge,len);

    for (int g=0; g<gpe.size(); g++){
        vertex mapped = GaussMapPointsEdge({gpe[g]}, edge);
        double uin  = use.eval(mapped, ml, "all", 
                      allwgts({gcellin[0],gcellin[1]}), gcellin, lu); 
        double uout = use.eval(mapped, ml, "all", 
                      allwgts({gcellout[0],gcellout[1]}), gcellout, lu); 

        //double LF = sqrt(vel.at(g)[0]*vel.at(g)[0] + vel.at(g)[1]*vel.at(g)[1]);
        double LF = abs(vel.at(g)[0]*unitNormal[0] + vel.at(g)[1]*unitNormal[1]);

        work += gwe[g] * LFflux(uin, uout, 
                                advfunc(uin,vel.at(g),unitNormal),
                                advfunc(uout,vel.at(g),unitNormal),LF) * len/2.0; 
    }

    return work;
}

double edgefluxintegral(const MeshInfo& mi, 
                        const indice& gcell,
                        const vertexSet& edge,
                        const Tensor<weights>& allwgts,
                        const vector<vertex>& vel,
                        multilevel& ml,
                        mluse& use,
                        double ** lu){

    double work = 0.0;

    //! Extract default gauess points and gauess weights.
    const valarray<double>& gwe = GaussWeightsEdge;
    const valarray<double>& gpe = GaussPointsEdge;

    // Get edge lendth and unit vector normal to the given edge
    double len = length(edge);
    vertex unitNormal = UnitNormal(edge,len);

    for (int g=0; g<gpe.size(); g++){
        vertex mapped = GaussMapPointsEdge({gpe[g]}, edge);
        double u  = use.eval(mapped, ml, "all", 
                    allwgts({gcell[0],gcell[1]}), gcell, lu); 

        //double LF = sqrt(vel.at(g)[0]*vel.at(g)[0] + vel.at(g)[1]*vel.at(g)[1]);
        double LF = abs(vel.at(g)[0]*unitNormal[0] + vel.at(g)[1]*unitNormal[1]);

        work += gwe[g] * LFflux(u, u, 
                                advfunc(u,vel.at(g),unitNormal),
                                advfunc(u,vel.at(g),unitNormal),LF) * len/2.0; 
    }

    return work;
}

double edgefluxintegral(const vertexSet& edge,
                        const vector<double>& bnval,
                        const vector<vertex>& vel){

    double work = 0.0;

    //! Extract default gauess points and gauess weights.
    const valarray<double>& gwe = GaussWeightsEdge;
    const valarray<double>& gpe = GaussPointsEdge;

    // Get edge lendth and unit vector normal to the given edge
    double len = length(edge);
    vertex unitNormal = UnitNormal(edge,len);

    for (int g=0; g<gpe.size(); g++){
        double val = bnval.at(g);

        //double LF = sqrt(vel.at(g)[0]*vel.at(g)[0] + vel.at(g)[1]*vel.at(g)[1]);
        double LF = abs(vel.at(g)[0]*unitNormal[0] + vel.at(g)[1]*unitNormal[1]);

        work += gwe[g] * LFflux(val, val, advfunc(val, vel.at(g), unitNormal),
                                          advfunc(val, vel.at(g), unitNormal),LF) *len/2.0;

    }

    return work;
}

double edgefluxintegral(const vertexSet& edge,
                        const double& bnval,
                        const vector<vertex>& vel){

    double work = 0.0;

    //! Extract default gauess points and gauess weights.
    const valarray<double>& gwe = GaussWeightsEdge;
    const valarray<double>& gpe = GaussPointsEdge;

    // Get edge lendth and unit vector normal to the given edge
    double len = length(edge);
    vertex unitNormal = UnitNormal(edge,len);

    for (int g=0; g<gpe.size(); g++){

        //double LF = sqrt(vel.at(g)[0]*vel.at(g)[0] + vel.at(g)[1]*vel.at(g)[1]);
        double LF = abs(vel.at(g)[0]*unitNormal[0] + vel.at(g)[1]*unitNormal[1]);

        work += gwe[g] * LFflux(bnval, bnval, advfunc(bnval, vel.at(g), unitNormal),
                                              advfunc(bnval, vel.at(g), unitNormal),LF) *len/2.0;

    }

    return work;
}

double getcellflux(const MeshInfo& mi, const indice& gcell,
                   const Tensor<double>& vertedge, 
                   const Tensor<double>& horiedge){

    double work = 0.0;

    double area = mi.cellArea.at(FlatIndic(mi,gcell));

    work += horiedge({gcell[0], gcell[1]});

    work -= horiedge({gcell[0], gcell[1]+1});

    work += vertedge({gcell[0], gcell[1]});

    work -= vertedge({gcell[0]+1, gcell[1]});

    work /= area;

    //cout << "At cell " << gcell[0] << "  " << gcell[1] << endl;
    //cout << "left   : " << horiedge({gcell[0], gcell[1]}) << " ";
    //cout << "right  : " << horiedge({gcell[0], gcell[1]+1}) << " ";
    //cout << "bottom : " << vertedge({gcell[0], gcell[1]}) << " ";
    //cout << "top    : " << vertedge({gcell[0]+1, gcell[1]}) << " ";
    //cout << endl << endl;;

    return work;  
}

// Similar way to obtain derivatives against u
int edgefluxintegral(const MeshInfo& mi, 
                     const indice& gcellin,
                     const indice& gcellout,
                     const vertexSet& edge,
                     const Tensor<weights>& allwgts,
                     const vector<vertex>& vel,
                     multilevel& ml,
                     mluse& use,
                     double ** lu,
                     derivative& der,
                     double& f){

    // Attention!!!!!! 
    // In serial code, stencil index equals global index
    // Which is not the case in parallel !!!!!!!!!!!
	 // fix it later

    f = 0.0;
    der.clear();

    //! Extract default gauess points and gauess weights.
    const valarray<double>& gwe = GaussWeightsEdge;
    const valarray<double>& gpe = GaussPointsEdge;

    // Get edge lendth and unit vector normal to the given edge
    double len = length(edge);
    vertex unitNormal = UnitNormal(edge,len);

    for (int g=0; g<gpe.size(); g++){

        vertex mapped = GaussMapPointsEdge({gpe[g]}, edge);

        derivative derin;
        use.der(mapped, ml, "all", gcellin, lu, mi, derin);

        //use.derpseudo(mapped, ml, "all", allwgts({gcellin[0], gcellin[1]}), 
        //              gcellin, lu, mi, derin);

        derivative derout;
        use.der(mapped, ml, "all", gcellout, lu, mi, derout);

        //use.derpseudo(mapped, ml, "all", allwgts({gcellout[0], gcellout[1]}), 
        //              gcellout, lu, mi, derout);


        double uin  = use.eval(mapped, ml, "all", 
                      allwgts({gcellin[0],gcellin[1]}), gcellin, lu); 
        double uout = use.eval(mapped, ml, "all", 
                      allwgts({gcellout[0],gcellout[1]}), gcellout, lu); 

        //double LF = sqrt(vel.at(g)[0]*vel.at(g)[0] + vel.at(g)[1]*vel.at(g)[1]);
        double LF = abs(vel.at(g)[0]*unitNormal[0] + vel.at(g)[1]*unitNormal[1]);

        f += gwe[g] * LFflux(uin, uout, 
                             advfunc(uin,vel.at(g),unitNormal),
                             advfunc(uout,vel.at(g),unitNormal),LF) * len/2.0; 
        derivative derfin;
        derivative derfout;
        derivative derLF;

        dadvfunc(derin , uin , vel.at(g), unitNormal, derfin);
        dadvfunc(derout, uout, vel.at(g), unitNormal, derfout);

//        cout << "derin : " << gcellin[0] << "  " << gcellin[1]<< endl;
//        unordered_map_print(derin);
//        cout << "derout : " << gcellout[0] << "  " << gcellout[1] << endl;
//        unordered_map_print(derout);
//        cout << "derfin : " << endl;
//        unordered_map_print(derfin);
//        cout << "derfout : " << endl;
//        unordered_map_print(derfout);

        derLFflux(derin, derout, derfin, derfout, LF, derLF);

        unordered_map_arithmetic(der, derLF, std::plus<double>(), 
                         gwe[g]*len/2.0, std::multiplies<double>());

    }

    return 1;
}

// Similar way to obtain derivatives against u
int edgefluxintegral(const MeshInfo& mi, 
                     const indice& gcell,
                     const vertexSet& edge,
                     const Tensor<weights>& allwgts,
                     const vector<vertex>& vel,
                     multilevel& ml,
                     mluse& use,
                     double ** lu,
                     derivative& der,
                     double& f){

    // Attention!!!!!! 
    // In serial code, stencil index equals global index
    // Which is not the case in parallel !!!!!!!!!!!
	 // fix it later

    f = 0.0;
    der.clear();

    //! Extract default gauess points and gauess weights.
    const valarray<double>& gwe = GaussWeightsEdge;
    const valarray<double>& gpe = GaussPointsEdge;

    // Get edge lendth and unit vector normal to the given edge
    double len = length(edge);
    vertex unitNormal = UnitNormal(edge,len);

    for (int g=0; g<gpe.size(); g++){

        vertex mapped = GaussMapPointsEdge({gpe[g]}, edge);

        derivative derin;
        use.der(mapped, ml, "all", gcell, lu, mi, derin);

        //use.derpseudo(mapped, ml, "all", allwgts({gcellin[0], gcellin[1]}), 
        //              gcellin, lu, mi, derin);

        double uin  = use.eval(mapped, ml, "all", 
                      allwgts({gcell[0],gcell[1]}), gcell, lu); 

        //double LF = sqrt(vel.at(g)[0]*vel.at(g)[0] + vel.at(g)[1]*vel.at(g)[1]);
        double LF = abs(vel.at(g)[0]*unitNormal[0] + vel.at(g)[1]*unitNormal[1]);

        f += gwe[g] * LFflux(uin, uin, 
                             advfunc(uin,vel.at(g),unitNormal),
                             advfunc(uin,vel.at(g),unitNormal),LF) * len/2.0; 
        derivative derfin;
        derivative derLF;

        dadvfunc(derin , uin , vel.at(g), unitNormal, derfin);

//        cout << "derin : " << gcellin[0] << "  " << gcellin[1]<< endl;
//        unordered_map_print(derin);
//        cout << "derout : " << gcellout[0] << "  " << gcellout[1] << endl;
//        unordered_map_print(derout);
//        cout << "derfin : " << endl;
//        unordered_map_print(derfin);
//        cout << "derfout : " << endl;
//        unordered_map_print(derfout);

        derLFflux(derin, derin, derfin, derfin, LF, derLF);

        unordered_map_arithmetic(der, derLF, std::plus<double>(), 
                         gwe[g]*len/2.0, std::multiplies<double>());
    }

    return 1;
}

int getcellflux(const MeshInfo& mi, const indice& gcell,
                const Tensor<double>& vertedge, 
                const Tensor<double>& horiedge,
                const Tensor<derivative>& vertedgeder,
                const Tensor<derivative>& horiedgeder,
                double& flux,
                derivative& dflux){

    double area = mi.cellArea.at(FlatIndic(mi,gcell));

    flux += horiedge({gcell[0], gcell[1]});

    flux -= horiedge({gcell[0], gcell[1]+1});

    flux += vertedge({gcell[0], gcell[1]});

    flux -= vertedge({gcell[0]+1, gcell[1]});

    flux /= area;

    // dflux 

    unordered_map_arithmetic(dflux, horiedgeder({gcell[0], gcell[1]}),
                             std::plus<double>());

    unordered_map_arithmetic(dflux, horiedgeder({gcell[0], gcell[1]+1}),
                             std::minus<double>());

    unordered_map_arithmetic(dflux, vertedgeder({gcell[0], gcell[1]}),
                             std::plus<double>());

    unordered_map_arithmetic(dflux, vertedgeder({gcell[0]+1, gcell[1]}),
                             std::minus<double>());

    unordered_map_arithmetic(dflux, 1.0/area,
                             std::multiplies<double>());

    return 1;
}
