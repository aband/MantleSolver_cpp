#ifndef RK_H_
#define RK_H_

#include "advectiveflux.h"
#include "transfunc.h"

int rk1(double dt, int Nt, Vec * insol, const MeshInfo& mi, 
        DM dmu, DM dmmesh,
        vector<reconstruction>& my_recon,
		  vector<tensorstencilpoly>& sten_lg,
		  vector<tensorstencilpoly>& sten_sm);

int rk2(double dt, int Nt, Vec * insol, const MeshInfo& mi, 
        DM dmu, DM dmmesh,
        vector<reconstruction>& my_recon,
		  vector<tensorstencilpoly>& sten_lg,
		  vector<tensorstencilpoly>& sten_sm);

#endif
