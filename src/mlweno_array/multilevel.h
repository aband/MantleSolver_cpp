#ifndef MULTILEVEL_H_
#define MULTILEVEL_H_

#include "tensorstencilpoly.h"

class singlelevel {

    public:

        singlelevel() {};

        ~singlelevel(){};

        int init(int order, int sizex, int sizey, const MeshInfo& mi);

        int myorder = 0;
        int mysizex = 0;
        int mysizey = 0;

        int totalsizex = 0;
        int totalsizey = 0;

        int totalsten = 0;

//    private:

        vector<tensorstencilpoly> stenpolys; 
};

class mlreconstruction {

    public:
        mlreconstruction() {};
        ~mlreconstruction() {};


        int prepare(const unordered_map<std::string, vector<indice>>& methods,
						  const unordered_map<std::string, indice>& size,
                    indice target, bool use_const, const MeshInfo& mi);

        int setWgts(double area, 
        const unordered_map<std::string,vector<double>>& allsigma);

        int printInfo();

    private:

        unordered_map<std::string, vector<indice>> my_method;

        unordered_map<std::string, vector<double>> linwgts;

        unordered_map<std::string, vector<double>> nonlinwgts;

};

class multilevel {

    public:

        multilevel() {};
        ~multilevel() {};

        int addLevel(int order, const std::string& name,
                     int sizex, int sizey,
                     const MeshInfo& mi);

        int updateSigma(double ** locvals);

        // Multilevel reconstruction 
        vector<mlreconstruction> my_recon;

        int printInfo();

    private:

        unordered_map<std::string, singlelevel> levels;
        unordered_map<std::string, vector<double>> allsigma;

};

#endif
