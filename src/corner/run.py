import myplot
import myanime_scalar
import myplot_scalar
import myplot_sampleedge
import myplot_single
import myplot_scalar_edge

#myplot.plot_velocity(20,20,1,'porochange','edgevel_stokes')
#myplot.plot_velocity(20,20,1,'porochange','edgevel_darcy')
#myplot_scalar.plot_scalar(40,40,1,'porochange','aveporo')
#myplot_scalar_edge.plot_edge_scalar(40,40,1,'porochange','porosity')
#myplot_scalar.plot_scalar(20,20,1,'porochange','newcellC')
#myplot_scalar.plot_scalar(20,20,1,'porochange','newcellH')
#myplot.plot_velocity(40,40,1,'porochange','edgevel_darcy_porosity')
#myplot.plot_velocity(40,40,1,'porochange','phaseave_vel')
#myplot_scalar_edge.plot_edge_scalar(20,20,1,'half','porosity')
#myplot.plot_velocity(20,20,1,'half','edgevel_darcy_porosity')
#myplot.plot_velocity(20,20,1,'half','phaseave_vel')

#myplot_scalar_edge.plot_edge_scalar(60,60,1,'porochange','porosity')
#myplot.plot_velocity(60,60,1,'porochange','edgevel_darcy_porosity')
#myplot.plot_velocity(60,60,1,'porochange','phaseave_vel')
#
#myplot_scalar.plot_scalar(20,20,1,'porochange','solidpotential')
#myplot_scalar.plot_scalar(20,20,1,'porochange','liquidpotential')

#myplot.plot_velocity(20,20,1,'half','edgevel_darcy_porosity')
#myplot.plot_velocity(20,20,1,'half','phaseave_vel')
#myplot_scalar.plot_scalar(20,20,1,'half','massconsv')

M=32
N=32
T=10
for i in range(T):
    k = i+1
    print(i)
    myplot_scalar.plot_scalar(M,N,k,'formalEvolve','massconsv')
    myplot_scalar_edge.plot_edge_scalar(M,N,k,'formalEvolve','porosity')
    myplot.plot_velocity(M,N,k,'formalEvolve','edgevel_darcy_porosity')
    myplot.plot_velocity(M,N,k,'formalEvolve','phaseave_vel')

#myplot_scalar.plot_scalar(20,20,1,'formalpreheat','cellC')
#myplot_scalar.plot_scalar(20,20,1,'formalpreheat','cellH')

#myplot_scalar.plot_scalar(20,20,1,'couplephase', 'newcellH')
#myplot_scalar.plot_scalar(20,20,1,'couplephase', 'newcellC')
#myplot_scalar.plot_scalar(20,20,1,'couplephase', 'aveporo')
#myplot_scalar.plot_scalar(20,20,1,'couplephase', 'massconsv')

#myplot.plot_velocity(40,40,1,'random','edgevel_darcy_porosity')
#myplot.plot_velocity(40,40,1,'random','phaseave_vel')
#myplot_scalar.plot_scalar(40,40,1,'random','aveporo')
#myplot_scalar.plot_scalar(40,40,1,'random','massconsv')

#myanime_scalar.anime_scalar(20,20,1,'porochange','cellC')

#myplot.plot_velocity(40,40,1,'half','edgevel_darcy_porosity')
#myplot.plot_velocity(40,40,1,'half','phaseave_vel')

#myplot.plot_velocity(50,50,1,'trans_alone','edgevel_stokes')

# Generate playable animation
#myanime_scalar.anime_scalar(25,25,1000,'trans_alone','cellH')
#myanime_scalar.anime_scalar(20,20,80,'trans_alone','cellC')

#myanime_scalar.anime_scalar(20,20,160,'advdiff','cellC')

#myanime_scalar.anime_scalar(50,50,1,'diff_alone','Init')
#myanime_scalar.anime_scalar(20,20,1000,'diff_alone','cellC')

#myplot_scalar.plot_scalar(5,5,1,'diff_alone','Init')

# Plot diffusion sampling points
#myplot_sampleedge.plot_scalar(4200, 'diff_alone')
#myplot_single.plot_single(400, 10,10,'diff_alone')

# Plot flux and cell averaged solutions
#myplot_scalar.plot_scalar(20,20,1,'diff_alone','cellC')
#myplot_scalar.plot_scalar(20,20,1,'diff_alone','flux')
#
#myplot_scalar.plot_scalar(20,20,999,'diff_alone','cellC')
#myplot_scalar.plot_scalar(20,20,999,'diff_alone','flux')

#myanime_scalar.anime_scalar(20,20,1000,'diff_alone','cellC')
#myanime_scalar.anime_scalar(50,50,80,'advdiff','cellC')

