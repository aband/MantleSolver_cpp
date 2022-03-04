

solution LaxFriedrichsFlux(const WenoMesh*& wm, double t, solution u_in, 
                           index_set& Stencillarge, vector<index_set>& StencilSmall point_index& target, 
                           vector<int>& Sorder, vector<int>& Lorder, point target_point, solution u_in){

    solution flux;
    solution u_out;

    switch (pos){
        case 0: // LEFT j,i-1
            u_out = WenoPointReconst(StencilLarge,StencilSmall,wm,{target[0]-1,target[1]},Sorder,Lorder,target_point); 
        case 1: // DOWN j-1,i
            u_out = WenoPointReconst(StencilLarge,StencilSmall,wm,{target[0],target[1]-1},Sorder,Lorder,target_point);
        case 2: // RIGHT j,i+1
            u_out = WenoPointReconst(StencilLarge,StencilSmall,wm,{target[0]+1,target[1]},Sorder,Lorder,target_point);
        case 3:
            u_out = WenoPointReconst(StencilLarge,StencilSmall,wm,{target[0],target[1]+1},Sorder,Lorder,target_point);
    }

    flux = (funcF(t,target_point,u_in) + funcF(t,target_point,u_out))*norm.x + (funcG(t,target_point,u_in) + funcG(t,target_point,u_out))*norm.y;

    // Determine stabilization parameter
    double alphaLF = MAX( pow(pow(dfuncF(t,target_point,u_in),2) + pow(dfuncF(t,target_point,u_out),2),0.5) , 
                          pow(pow(dfuncG(t,target_point,u_in),2) + pow(dfuncG(t,target_point,u_out),2),0.5) );


    return flux;
}
