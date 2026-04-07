import myplot
import myanime_scalar
import myplot_scalar
import myplot_sampleedge
import myplot_single

#myplot.plot_velocity(20,20,1,'advdiff','edgevel_stokes')
#myplot.plot_velocity(40,40,1,'half','edgevel_darcy')
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
myplot_scalar.plot_scalar(20,20,1,'diff_alone','cellC')
myplot_scalar.plot_scalar(20,20,1,'diff_alone','flux')

myplot_scalar.plot_scalar(20,20,999,'diff_alone','cellC')
myplot_scalar.plot_scalar(20,20,999,'diff_alone','flux')

#myanime_scalar.anime_scalar(20,20,1000,'diff_alone','cellC')
myanime_scalar.anime_scalar(20,20,16,'advdiff','cellC')

