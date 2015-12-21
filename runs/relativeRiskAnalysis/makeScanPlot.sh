#!/bin/bash

###############################################################################
#
#    Particle Therapy Project Bergen (PTPB) - tools and models for research in
#    cancer therapy using particle beams.
#
#    Copyright (C) 2015 Particle Therapy Group Bergen
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
###############################################################################

if test -n "$1" ; then
    TECH_ORGAN_NAME="$1"
else
    TECH_ORGAN_NAME="cion_bladder"
fi

if echo "$TECH_ORGAN_NAME" | grep '^cion_' > /dev/null ; then
    TITLE_DENOM="C-ion"
else
    TITLE_DENOM="IMPT"
fi

create_datafile () {
octave -q <<EOF
z = $1;
w = $2;
load('$3');
datafile = '$4';
errmul = $5;
[fid, msg] = fopen(datafile, 'w');
if fid == -1
    error('Could not open output file %s: %s', datafile, msg);
end
[Nx, Ny, Nz, Nw] = size(RR);
for x = 1:Nx
    for y = 1:Ny
        X = alpha(x, y, z, w);
        Y = beta(x, y, z, w);
        Z = RR(x, y, z, w) + errmul .* RRerr(x, y, z, w);
        fprintf(fid, '%g\t%g\t%g\n', X, Y, Z);
    end
    fprintf(fid, '\n');
end
EOF
}

get_rbe_min () {
octave -q <<EOF
z = $1;
w = $2;
load('$3');
printf('%g', RBEmin(1,1,z,w));
EOF
}

get_rbe_max () {
octave -q <<EOF
z = $1;
w = $2;
load('$3');
printf('%g', RBEmax(1,1,z,w));
EOF
}

HYPERCUBEFILE="hypercube_${TECH_ORGAN_NAME}.mat"

FILENUMS=""
for N in 1 2 3 4; do
    for M in 1 2 3 4; do
        FILENUMS="$FILENUMS $N$M"
        create_datafile $N $M "$HYPERCUBEFILE" "relative_risk_${TECH_ORGAN_NAME}_$N$M.dat" 0
        create_datafile $N $M "$HYPERCUBEFILE" "relative_risk_high_${TECH_ORGAN_NAME}_$N$M.dat" +1
        create_datafile $N $M "$HYPERCUBEFILE" "relative_risk_low_${TECH_ORGAN_NAME}_$N$M.dat" -1
    done
done

for N in $FILENUMS ; do
    gnuplot <<EOF
set contour base
set cntrparam levels incr 0,0.25,3
set cntrparam cubicspline
set cntrparam order 7
set cntrparam points 10
set xrange [0.01:0.79]
set yrange [0.01:0.05]
set zrange [0.25:5]
set cbrange [0.25:2.25]
unset surface
set table "contours_${TECH_ORGAN_NAME}_$N.dat"
splot "relative_risk_${TECH_ORGAN_NAME}_$N.dat"
set table "contours_high_${TECH_ORGAN_NAME}_$N.dat"
splot "relative_risk_high_${TECH_ORGAN_NAME}_$N.dat"
set table "contours_low_${TECH_ORGAN_NAME}_$N.dat"
splot "relative_risk_low_${TECH_ORGAN_NAME}_$N.dat"
EOF

    ./makeContourLabels.py contours_${TECH_ORGAN_NAME}_$N.dat \
                           processed_contours_${TECH_ORGAN_NAME}_$N.dat \
                           contour_labels_${TECH_ORGAN_NAME}_$N.dat
done

RBEMIN_1=`get_rbe_min 1 1 "$HYPERCUBEFILE"`
RBEMIN_2=`get_rbe_min 2 1 "$HYPERCUBEFILE"`
RBEMIN_3=`get_rbe_min 3 1 "$HYPERCUBEFILE"`
RBEMIN_4=`get_rbe_min 4 1 "$HYPERCUBEFILE"`
RBEMAX_1=`get_rbe_max 1 1 "$HYPERCUBEFILE"`
RBEMAX_2=`get_rbe_max 1 2 "$HYPERCUBEFILE"`
RBEMAX_3=`get_rbe_max 1 3 "$HYPERCUBEFILE"`
RBEMAX_4=`get_rbe_max 1 4 "$HYPERCUBEFILE"`

gnuplot <<EOF
set terminal postscript eps color enhanced dashed size 16cm,16cm
set output "${TECH_ORGAN_NAME}.eps"

set multiplot layout 4,4 columnsfirst downwards scale 1, 1

unset key
set pm3d map
set pm3d interpolate 4,4
# set title "Cion - Bladder - RBEmin = 1.05, RBEmax = 2"
set xrange [0.01:0.79]
set yrange [0.01:0.05]
set zrange [0.25:5]
set cbrange [0.25:2.25]

px1 = 0.085
px5 = 0.965
px2 = (px5 - px1) / 4 * 1 + px1
px3 = (px5 - px1) / 4 * 2 + px1
px4 = (px5 - px1) / 4 * 3 + px1

py1 = 0.14
py5 = 0.965
py2 = (py5 - py1) / 4 * 1 + py1
py3 = (py5 - py1) / 4 * 2 + py1
py4 = (py5 - py1) / 4 * 3 + py1

xlabelvalue = "{/Symbol a}  [Gy^{-1}]"
ylabelvalue = "{/Symbol b}  [Gy^{-2}]"

###############################################################################

set title "RBE_{min} = $RBEMIN_1"
unset colorbox
unset xlabel
set format x ''
unset cblabel
set ylabel ylabelvalue
set ytics 0.01,0.005,0.05
set lmargin at screen px1
set rmargin at screen px2
set bmargin at screen py4
set tmargin at screen py5
splot "relative_risk_${TECH_ORGAN_NAME}_11.dat", \
    "processed_contours_${TECH_ORGAN_NAME}_11.dat" using 1:2:3:3 with lines lw 2 lc rgb "#000000", \
    "contour_labels_${TECH_ORGAN_NAME}_11.dat" using 1:2:3:3 with labels

unset title
set ytics 0.01,0.005,0.045
set lmargin at screen px1
set rmargin at screen px2
set bmargin at screen py3
set tmargin at screen py4
splot "relative_risk_${TECH_ORGAN_NAME}_12.dat", \
    "processed_contours_${TECH_ORGAN_NAME}_12.dat" using 1:2:3:3 with lines lw 2 lc rgb "#000000", \
    "contour_labels_${TECH_ORGAN_NAME}_12.dat" using 1:2:3:3 with labels

set lmargin at screen px1
set rmargin at screen px2
set bmargin at screen py2
set tmargin at screen py3
splot "relative_risk_${TECH_ORGAN_NAME}_13.dat", \
    "processed_contours_${TECH_ORGAN_NAME}_13.dat" using 1:2:3:3 with lines lw 2 lc rgb "#000000", \
    "contour_labels_${TECH_ORGAN_NAME}_13.dat" using 1:2:3:3 with labels

set xlabel xlabelvalue
set format x
set lmargin at screen px1
set rmargin at screen px2
set bmargin at screen py1
set tmargin at screen py2
splot "relative_risk_${TECH_ORGAN_NAME}_14.dat", \
    "processed_contours_${TECH_ORGAN_NAME}_14.dat" using 1:2:3:3 with lines lw 2 lc rgb "#000000", \
    "contour_labels_${TECH_ORGAN_NAME}_14.dat" using 1:2:3:3 with labels

###############################################################################

set title "RBE_{min} = $RBEMIN_2"
unset xlabel
set format x ''
unset ylabel
set format y ''
set lmargin at screen px2
set rmargin at screen px3
set bmargin at screen py4
set tmargin at screen py5
splot "relative_risk_${TECH_ORGAN_NAME}_21.dat", \
    "processed_contours_${TECH_ORGAN_NAME}_21.dat" using 1:2:3:3 with lines lw 2 lc rgb "#000000", \
    "contour_labels_${TECH_ORGAN_NAME}_21.dat" using 1:2:3:3 with labels

unset title
set lmargin at screen px2
set rmargin at screen px3
set bmargin at screen py3
set tmargin at screen py4
splot "relative_risk_${TECH_ORGAN_NAME}_22.dat", \
    "processed_contours_${TECH_ORGAN_NAME}_22.dat" using 1:2:3:3 with lines lw 2 lc rgb "#000000", \
    "contour_labels_${TECH_ORGAN_NAME}_22.dat" using 1:2:3:3 with labels

set lmargin at screen px2
set rmargin at screen px3
set bmargin at screen py2
set tmargin at screen py3
splot "relative_risk_${TECH_ORGAN_NAME}_23.dat", \
    "processed_contours_${TECH_ORGAN_NAME}_23.dat" using 1:2:3:3 with lines lw 2 lc rgb "#000000", \
    "contour_labels_${TECH_ORGAN_NAME}_23.dat" using 1:2:3:3 with labels

set xlabel xlabelvalue
set format x
set lmargin at screen px2
set rmargin at screen px3
set bmargin at screen py1
set tmargin at screen py2
splot "relative_risk_${TECH_ORGAN_NAME}_24.dat", \
    "processed_contours_${TECH_ORGAN_NAME}_24.dat" using 1:2:3:3 with lines lw 2 lc rgb "#000000", \
    "contour_labels_${TECH_ORGAN_NAME}_24.dat" using 1:2:3:3 with labels

###############################################################################

set title "RBE_{min} = $RBEMIN_3"
unset xlabel
set format x ''
unset ylabel
set format y ''
set lmargin at screen px3
set rmargin at screen px4
set bmargin at screen py4
set tmargin at screen py5
splot "relative_risk_${TECH_ORGAN_NAME}_31.dat", \
    "processed_contours_${TECH_ORGAN_NAME}_31.dat" using 1:2:3:3 with lines lw 2 lc rgb "#000000", \
    "contour_labels_${TECH_ORGAN_NAME}_31.dat" using 1:2:3:3 with labels

unset title
set lmargin at screen px3
set rmargin at screen px4
set bmargin at screen py3
set tmargin at screen py4
splot "relative_risk_${TECH_ORGAN_NAME}_32.dat", \
    "processed_contours_${TECH_ORGAN_NAME}_32.dat" using 1:2:3:3 with lines lw 2 lc rgb "#000000", \
    "contour_labels_${TECH_ORGAN_NAME}_32.dat" using 1:2:3:3 with labels

set lmargin at screen px3
set rmargin at screen px4
set bmargin at screen py2
set tmargin at screen py3
splot "relative_risk_${TECH_ORGAN_NAME}_33.dat", \
    "processed_contours_${TECH_ORGAN_NAME}_33.dat" using 1:2:3:3 with lines lw 2 lc rgb "#000000", \
    "contour_labels_${TECH_ORGAN_NAME}_33.dat" using 1:2:3:3 with labels

set xlabel xlabelvalue
set format x
set lmargin at screen px3
set rmargin at screen px4
set bmargin at screen py1
set tmargin at screen py2
splot "relative_risk_${TECH_ORGAN_NAME}_34.dat", \
    "processed_contours_${TECH_ORGAN_NAME}_34.dat" using 1:2:3:3 with lines lw 2 lc rgb "#000000", \
    "contour_labels_${TECH_ORGAN_NAME}_34.dat" using 1:2:3:3 with labels

###############################################################################

pRBElabelx = px5 + 0.02
pRBElabely1 = 0.5*(py5 + py4)
pRBElabely2 = 0.5*(py4 + py3)
pRBElabely3 = 0.5*(py3 + py2)
pRBElabely4 = 0.5*(py2 + py1)

set title "RBE_{min} = $RBEMIN_4"
unset xlabel
set format x ''
unset ylabel
set format y ''
set label 1 "RBE_{max} = $RBEMAX_1" at screen pRBElabelx, pRBElabely1 center rotate
set lmargin at screen px4
set rmargin at screen px5
set bmargin at screen py4
set tmargin at screen py5
splot "relative_risk_${TECH_ORGAN_NAME}_41.dat", \
    "processed_contours_${TECH_ORGAN_NAME}_41.dat" using 1:2:3:3 with lines lw 2 lc rgb "#000000", \
    "contour_labels_${TECH_ORGAN_NAME}_41.dat" using 1:2:3:3 with labels

unset title
set label 1 "RBE_{max} = $RBEMAX_2" at screen pRBElabelx, pRBElabely2 center rotate
set lmargin at screen px4
set rmargin at screen px5
set bmargin at screen py3
set tmargin at screen py4
splot "relative_risk_${TECH_ORGAN_NAME}_42.dat", \
    "processed_contours_${TECH_ORGAN_NAME}_42.dat" using 1:2:3:3 with lines lw 2 lc rgb "#000000", \
    "contour_labels_${TECH_ORGAN_NAME}_42.dat" using 1:2:3:3 with labels

set label 1 "RBE_{max} = $RBEMAX_3" at screen pRBElabelx, pRBElabely3 center rotate
set lmargin at screen px4
set rmargin at screen px5
set bmargin at screen py2
set tmargin at screen py3
splot "relative_risk_${TECH_ORGAN_NAME}_43.dat", \
    "processed_contours_${TECH_ORGAN_NAME}_43.dat" using 1:2:3:3 with lines lw 2 lc rgb "#000000", \
    "contour_labels_${TECH_ORGAN_NAME}_43.dat" using 1:2:3:3 with labels

cbwidth = px5 - px1
set colorbox horizontal user origin px1,0.055 size cbwidth,0.013
set cblabel "Relative Risk (VMAT/${TITLE_DENOM})"
set xlabel xlabelvalue
set format x
set label 1 "RBE_{max} = $RBEMAX_4" at screen pRBElabelx, pRBElabely4 center rotate
set lmargin at screen px4
set rmargin at screen px5
set bmargin at screen py1
set tmargin at screen py2
splot "relative_risk_${TECH_ORGAN_NAME}_44.dat", \
    "processed_contours_${TECH_ORGAN_NAME}_44.dat" using 1:2:3:3 with lines lw 2 lc rgb "#000000", \
    "contour_labels_${TECH_ORGAN_NAME}_44.dat" using 1:2:3:3 with labels
EOF
