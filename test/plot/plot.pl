set title "Burger's equation with Oblique Initial Condition"
#plot 'build/Pressure.data' matrix with image notitle

#plot 'Pressure.data' matrix with image notitle

reset
set contour
unset surface
set cntrparam level incremental -1.0, 0.5, 1.0

set view map

set dgrid3d 100,100,4

set table 'contour.data'
splot 'Pressure.data' matrix
unset table

unset contour
set surface 
set table 'grid.data'
splot 'Pressure.data' matrix
unset table

reset
set pm3d map
unset key

set palette defined (0 '#352a87', 1 '#0363e1',2 '#1485d4', 3 '#06a7c6', 4 '#38b99e', 5 '#92bf73', 6 '#d9ba56', 7 '#fcce2e', 8 '#f9fb0e')
set autoscale fix
set grid

splot 'grid.data' w pm3d, 'contour.data' w l lc rgb 'black'

pause -1
