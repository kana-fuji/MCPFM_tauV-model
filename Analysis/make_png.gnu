#run
#gnuplot -e "No='$1';tmax=$2" pf_cell_2D.gnu

stats "../DATA/".No."/timestamp.dat" using 1 nooutput
tmax=int(STATS_max)
#tmax=420000
#print tmax

#stats "../IN/".No.".param" using 2 every ::2::2 nooutput
dt=1.0
#print dt

#stats "../IN/".No.".param" using 3 every ::2::2 nooutput
#samplingdt=STATS_max
#print samplingdt

ttmin=0
tmin=0
tstep=10

#tstep=1.0/dt*samplingdt
#tstep=1000
#print tstep

imin=0
imax=5.0

#gif
#set encoding utf8
#set terminal gif animate delay 6 optimize size 640,480
#set output "DATA/GIF/pf_cell_2D_".No.".gif"

set pm3d map
set palette maxcolors 50

unset colorbox

set size square
set xrange [0:1024]
set yrange [0:1024]
#set yrange [] reverse
set zrange [-100.0:10.0]
#set cbrange [0.1:0.3]
#set cbrange [-1.0:1.0]
#set xlabel "x"
#set ylabel "y"
#set zlabel "u"
#set xtics ("0" 0, "2" 100, "4" 200, "6" 300, "8" 400, "10" 500, "12" 600, "14" 700)
#set ytics ("0" 767, "2" 668, "4" 568, "6" 468, "8" 368, "10" 268,"12" 168, "14" 68)
#set ytics ("0" 512, "2" 412, "4" 312, "6" 212, "8" 112, "10" 12)

set noxtics
set noytics

#set palette rgbformulae 22,13,-31
#set palette defined(0.0"white",0.5"black",1.0"red")
#set palette maxcolors 3
#set palette rgbformulae 22,13,-31
#set palette defined(0"white",0.2"black",1"red")

#set palette defined(-1"blue",-0.51"blue",-0.5"black",-0.49"white",0.0"white",0.19"white",0.2"black",0.21"black",1"red")

#do for [t=tmin:tmax:tstep] {
#   time=sprintf("%.2f",t*dt) 
#   set title "t=".time
#   sp "DATA/".No."/u_".t.".dat" matrix not w pm3d
#}
#unset output

#-----------------

set pm3d map
set palette maxcolors 50
set terminal png size 480,480
set output "PNG/".No.".png"	

set cbrange [0.0:1.0]
set lmargin screen 0.05
set rmargin screen 0.95
set tmargin screen 0.95
set bmargin screen 0.05


t=tmax

time=sprintf("%.2f",t*dt) 
set title "t=".time

set palette defined(0.0"white",0.5"black",1.0"red")
sp "../DATA/".No."/u_".t.".dat" matrix not w pm3d

#set multiplot layout 1,2
#    set palette defined(0.0"white",0.5"black",1.0"red")
#    sp "../DATA/".No."/u_".t.".dat" matrix not w pm3d

#    set palette defined(0.0"white",0.5"black",1.0"blue")
#    sp "../DATA/".No."/s_".t.".dat" matrix not w pm3d
#unset multiplot

#    sp "../DATA/".No."/c_".t.".dat" matrix not w pm3d
#    set palette gray
#    sp "../DATA/".No."/p_".t.".dat" matrix not w pm3d

#-----------------

#set terminal gif animate delay 6 optimize size 640,480
#set output "DATA/GIF/ecm_2D_".No.".gif"

#set cbrange [0.0:1.0]
#set palette defined(0"white",0.5"black",1"yellow")
#set palette defined(0"white",0.48"white",0.5"black",0.52"yellow",1"yellow")

#do for [t=tmin:tmax:tstep] {
#   time=sprintf("%.2f",t*dt) 
#   set title "t=".time
#   sp "DATA/".No."/c_".t.".dat" matrix not w pm3d
#}
#unset output

#-----------------


#set terminal gif animate delay 3 optimize size 1280,480
#set output "GIF/p_2D_".No.".gif"
#set output "DATA/GIF/uscp_2D_".No.".gif"


#set cbrange [0.0:1.0]
#set palette defined(0"white",0.5"black",1"yellow")
#set palette defined(0"white",0.48"white",0.5"black",0.52"yellow",1"yellow")


#do for [t=tmin:tmax:tstep] {
#   time=sprintf("%.2f",t*dt) 
#   set title "t=".time
#   set palette gray
#   sp "../DATA/".No."/p_".t.".dat" matrix not w pm3d
#}
#unset output

#-----------------


#set terminal gif animate delay 2 optimize size 1280,480
#set output "GIF/uscp_2D_".No.".gif"

#set cbrange [0.0:1.0]
#set palette defined(0"white",0.5"black",1"yellow")
#set palette defined(0"white",0.48"white",0.5"black",0.52"yellow",1"yellow")

#do for [t=tmin:tmax:tstep] {

#t=tmin;
#while(t<tmax){

#   time=sprintf("%.2f",t*dt) 
#   set title "t=".time
#   set multiplot layout 1,4
      
#set palette defined(0.0"white",0.19"white",0.2"black",0.21"black",1"red")

#      set palette rgbformulae 22,13,-31
#      set zrange [-100.0:10.0]
#      sp "../DATA/".No."/u_".t.".dat" matrix not w pm3d
#      set zrange [0.00001:10.0]
#      sp "../DATA/".No."/s_".t.".dat" matrix not w pm3d
#      set zrange [-100.0:10.0]
#      sp "../DATA/".No."/c_".t.".dat" matrix not w pm3d
#      set palette gray
#      sp "../DATA/".No."/p_".t.".dat" matrix not w pm3d

#   unset multiplot
#   t=t+tstep
#}
#unset output


