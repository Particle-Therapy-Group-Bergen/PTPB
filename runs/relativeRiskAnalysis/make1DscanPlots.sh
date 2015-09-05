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
    DATAFILE="$1"
else
    DATAFILE="alpha_scan_cion_bladder_sample_data.mat"
fi

if test -n "$2" ; then
    OUTPUTFILE="$2"
else
    OUTPUTFILE="alpha_scan_cion_bladder.dat"
fi

if test -n "$3" ; then
    BINS="$3"
else
    BINS="0.01:0.02:0.8"
fi

if test -n "$4" ; then
    DIM="$4"
else
    DIM=2
fi

if test -n "$5" ; then
    PLOTFILE="$5"
else
    PLOTFILE="alpha_scan_cion_bladder.eps"
fi

if test -n "$6" ; then
    XLABEL="$6"
else
    XLABEL="{/Symbol a}  [Gy^{-1}]"
fi


octave -q <<EOF
datafile = '$DATAFILE';
outputfile = '$OUTPUTFILE';
load(datafile);
[fid, msg] = fopen(outputfile, 'w');
if fid == -1
    error('Could not open output file %s: %s', datafile, msg);
end
bins = $BINS;
for n = 1:length(bins)-1
    low = bins(n);
    high = bins(n+1);
    index = find(low <= Results(:,$DIM) & Results(:,$DIM) < high);
    data = Results(index,1);
    Q = quantile(data, [0.025, 0.5, 0.975], 1, 8);
    M = mean(data);
    x = (low + high) * 0.5;
    fprintf(fid, '%g\t%g\t%g\t%g\t%g\t%d\n', x, M, Q(1), Q(2), Q(3), length(data));
end
EOF


gnuplot <<EOF
set terminal postscript eps color enhanced dashed size 8cm,6cm
set output "$PLOTFILE"
unset key
set xlabel "$XLABEL"
set ylabel "Relative Risk"
plot "$OUTPUTFILE" using 1:2:3:5 with errorbars
EOF
