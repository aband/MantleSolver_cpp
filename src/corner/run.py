import myplot
import myanime_scalar
import myplot_scalar

#myplot.plot_velocity(40,40,1,'half','edgevel_stokes')
#myplot.plot_velocity(40,40,1,'half','edgevel_darcy')
#myplot.plot_velocity(40,40,1,'half','edgevel_darcy_porosity')
#myplot.plot_velocity(40,40,1,'half','phaseave_vel')

#myplot.plot_velocity(50,50,1,'trans_alone','edgevel_stokes')

# Generate playable animation
#myanime_scalar.anime_scalar(25,25,1000,'trans_alone','cellH')
myanime_scalar.anime_scalar(20,20,80,'trans_alone','cellC')

#myplot_scalar.plot_scalar(5,5,1,'trans_alone','fluxH')
