# gnuplot script file to plot scattered temperature distribution
# over the 2D domain

reset
#set size square
set key off
set xrange [0:2]
set yrange [0:1]
set format x ""
set format y ""
set palette rgbformulae 33,13,10
#set palette maxcolors 11
set cbrange [206:410]
set view map
splot 'test_case_4.dat' with points palette
#set dgrid3d 100,200
#set hidden3d
#set pm3d map
#splot 'data.dat'
