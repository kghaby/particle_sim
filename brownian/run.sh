#!/bin/bash
set -e

###what file am i writing
f="traj.dat"
### delete files if they exist
rm $f


### simulate 
#python3 atom.py >> traj.dat
#python3 atom_x.py >> traj_x.dat
python3 brownian.py >> traj.dat

echo "successfully simulated woohoo"

#sed -i '' -e '$ d' traj.dat
### plot

range_x=$(bc<<<"scale=5; $(awk 'BEGIN{x=   0}{if ($1>0+x) x=$1} END{print x}' $f)- $(awk 'BEGIN{x=   0}{if ($1<0+x) x=$1} END{print x}' $f)" )
range_y=$(bc<<<"scale=5; $(awk 'BEGIN{y=   0}{if ($2>0+y) y=$2} END{print y}' $f)- $(awk 'BEGIN{y=   0}{if ($2<0+y) y=$2} END{print y}' $f)" )
range_z=$(bc<<<"scale=5; $(awk 'BEGIN{z=   0}{if ($3>0+z) z=$3} END{print z}' $f)- $(awk 'BEGIN{z=   0}{if ($3<0+z) z=$3} END{print z}' $f)" )
	#max value minus min value
bw_x=$(bc<<<"scale=5; ${range_x#-}/100")
bw_y=$(bc<<<"scale=5; ${range_y#-}/100")
bw_z=$(bc<<<"scale=5; ${range_z#-}/100")
	#divide absolute value by 100
echo "range_x was $range_x, so binwidth_x was $bw_x"
echo "range_y was $range_y, so binwidth_y was $bw_y"
echo "range_z was $range_z, so binwidth_z was $bw_z"

gnuplot -persist << EOF

set key t r

set term x11 1
set xrange [-20:20]
f(x)=1*(0.0001*x**4 - 0.005*x**2 - 2*exp(-1*(x+5)**2) -3*exp(-1*(x-5)**2))
g(x)=-1*(0.0004*x**3 - 0.01*x + (4 *(x+5))*exp(-1*(x+5)**2) + (6*(x-5))*exp(-1*(x-5)**2))
set title "force"
set xzeroaxis ls -1
plot f(x) w l title "potential energy",\
g(x) w l title "force"

set xrange [*:*]
set yrange [*:*]
set zrange [*:*]


set term x11 2 
#set output "/tmp/plot1.png"
set title "position over time for each dimension"
set xlabel "Time"
set ylabel "Position"
plot "$f" using 1 with lines title "x",\
"$f" using 2 with lines title "y",\
"$f" using 3 with lines title "z"

set term x11 3
#set output "/tmp/plot2.png"
set title "histogram"
set style fill transparent
binwidth_x=$bw_x
binwidth_y=$bw_y
binwidth_z=$bw_z
bin(x,width)=width*floor(x/width)
plot "$f" u (bin(\$1,binwidth_x)):(1.0) smooth frequency w boxes title "x",\
"$f" u (bin(\$2,binwidth_y)):(1.0) smooth frequency w boxes title "y",\
"$f" u (bin(\$3,binwidth_z)):(1.0) smooth frequency w boxes title "z"

set term x11 4
set title "3d motion"
splot "$f" u 1:2:3 w l

EOF


python3 animate.py


