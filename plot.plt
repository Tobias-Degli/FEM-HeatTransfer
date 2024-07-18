set term wxt enhanced persist
set xlabel 'X'
set ylabel 'Z'
set autoscale xy

set pm3d map interpolate 20,20 
set pal maxcolors 7
#set palette rbgformulae
#set palette gray #colorgray
splot 	'dados.txt' using 1:2:3 notitle

