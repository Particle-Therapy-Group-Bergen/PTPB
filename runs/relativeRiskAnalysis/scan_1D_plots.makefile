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

# This generates the scan plots after the *_scan_samples.makefile files have been
# used to create the *.mat input files.

TECH_ORGAN_TYPES = cion_bladder cion_rectum proton_bladder proton_rectum

PLOT_FILE_TMPL = alpha_scan_TECHORGAN.eps \
                 beta_scan_TECHORGAN.eps \
                 RBEmin_scan_TECHORGAN.eps \
                 RBEmax_scan_TECHORGAN.eps

PLOT_FILES = $(foreach N,$(TECH_ORGAN_TYPES),$(subst TECHORGAN,$(N),$(PLOT_FILE_TMPL)))

OUTPUT_FILES = $(subst .eps,.dat,$(PLOT_FILES))

.PHONY: all clean cleanall

all: $(PLOT_FILES) summary_1d_plots.eps summary_cion_1d_plots.eps summary_proton_1d_plots.eps

clean:
	rm -rf $(OUTPUT_FILES)

cleanall: clean
	rm -rf $(PLOT_FILES)


alpha_scan_%.eps: alpha_scan_%_sample_data.mat scan_1D_plots.makefile make1DscanPlots.sh
	./make1DscanPlots.sh $< $(subst .eps,.dat,$@) '0.01:0.02:0.8' 2 $@ '{/Symbol a}  [Gy^{-1}]'

beta_scan_%.eps: beta_scan_%_sample_data.mat scan_1D_plots.makefile make1DscanPlots.sh
	./make1DscanPlots.sh $< $(subst .eps,.dat,$@) '0.01:0.0025:0.05' 3 $@ '{/Symbol b}  [Gy^{-2}]'

RBEmin_scan_cion_%.eps: RBEmin_scan_cion_%_sample_data.mat scan_1D_plots.makefile make1DscanPlots.sh
	./make1DscanPlots.sh $< $(subst .eps,.dat,$@) '1.0:0.01:1.4' 4 $@ 'RBE_{min}  [Gy]'

RBEmin_scan_proton_%.eps: RBEmin_scan_proton_%_sample_data.mat scan_1D_plots.makefile make1DscanPlots.sh
	./make1DscanPlots.sh $< $(subst .eps,.dat,$@) '1.0:0.005:1.1' 4 $@ 'RBE_{min}  [Gy]'

RBEmax_scan_cion_%.eps: RBEmax_scan_cion_%_sample_data.mat scan_1D_plots.makefile make1DscanPlots.sh
	./make1DscanPlots.sh $< $(subst .eps,.dat,$@) '1:0.2:9' 5 $@ 'RBE_{max}  [Gy]'

RBEmax_scan_proton_%.eps: RBEmax_scan_proton_%_sample_data.mat scan_1D_plots.makefile make1DscanPlots.sh
	./make1DscanPlots.sh $< $(subst .eps,.dat,$@) '1.05:0.01:1.6' 5 $@ 'RBE_{max}  [Gy]'


define gnuplot_1d_summary_plot_script
set terminal postscript eps color enhanced dashed size 16cm,24cm
set output "summary_1d_plots.eps"
set multiplot layout 4,2 columnsfirst downwards
set ylabel "Relative Risk"
set xlabel "{/Symbol a}  [Gy^{-1}]"
set xrange [0:0.8]
set yrange [0:3.2]
set key at 0.75,3
plot "alpha_scan_cion_bladder.dat" using 1:2 with lines lt 1 lw 2 lc rgb "#CE0000" title "VMAT/C-ion bladder", \
     "alpha_scan_cion_bladder.dat" using 1:3 with lines lt 3 lw 1 lc rgb "#CE0000" notitle, \
     "alpha_scan_cion_bladder.dat" using 1:5 with lines lt 3 lw 1 lc rgb "#CE0000" notitle, \
     "alpha_scan_cion_rectum.dat" using 1:2 with lines lt 1 lw 2 lc rgb "#1080FF" title "VMAT/C-ion rectum", \
     "alpha_scan_cion_rectum.dat" using 1:3 with lines lt 3 lw 1 lc rgb "#1080FF" notitle, \
     "alpha_scan_cion_rectum.dat" using 1:5 with lines lt 3 lw 1 lc rgb "#1080FF" notitle
unset key
set xlabel "{/Symbol b}  [Gy^{-2}]"
set xrange [0.01:0.05]
set yrange [0:3.2]
plot "beta_scan_cion_bladder.dat" using 1:2 with lines lt 1 lw 2 lc rgb "#CE0000" title "VMAT/C-ion bladder", \
     "beta_scan_cion_bladder.dat" using 1:3 with lines lt 3 lw 1 lc rgb "#CE0000" notitle, \
     "beta_scan_cion_bladder.dat" using 1:5 with lines lt 3 lw 1 lc rgb "#CE0000" notitle, \
     "beta_scan_cion_rectum.dat" using 1:2 with lines lt 1 lw 2 lc rgb "#1080FF" title "VMAT/C-ion rectum", \
     "beta_scan_cion_rectum.dat" using 1:3 with lines lt 3 lw 1 lc rgb "#1080FF" notitle, \
     "beta_scan_cion_rectum.dat" using 1:5 with lines lt 3 lw 1 lc rgb "#1080FF" notitle
set xlabel "RBE_{max}"
set xrange [1:9]
set yrange [0:3.2]
plot "RBEmax_scan_cion_bladder.dat" using 1:2 with lines lt 1 lw 2 lc rgb "#CE0000" title "VMAT/C-ion bladder", \
     "RBEmax_scan_cion_bladder.dat" using 1:3 with lines lt 3 lw 1 lc rgb "#CE0000" notitle, \
     "RBEmax_scan_cion_bladder.dat" using 1:5 with lines lt 3 lw 1 lc rgb "#CE0000" notitle, \
     "RBEmax_scan_cion_rectum.dat" using 1:2 with lines lt 1 lw 2 lc rgb "#1080FF" title "VMAT/C-ion rectum", \
     "RBEmax_scan_cion_rectum.dat" using 1:3 with lines lt 3 lw 1 lc rgb "#1080FF" notitle, \
     "RBEmax_scan_cion_rectum.dat" using 1:5 with lines lt 3 lw 1 lc rgb "#1080FF" notitle
set xlabel "RBE_{min}"
set xrange [1:1.4]
set yrange [0:3.2]
plot "RBEmin_scan_cion_bladder.dat" using 1:2 with lines lt 1 lw 2 lc rgb "#CE0000" title "VMAT/C-ion bladder", \
     "RBEmin_scan_cion_bladder.dat" using 1:3 with lines lt 3 lw 1 lc rgb "#CE0000" notitle, \
     "RBEmin_scan_cion_bladder.dat" using 1:5 with lines lt 3 lw 1 lc rgb "#CE0000" notitle, \
     "RBEmin_scan_cion_rectum.dat" using 1:2 with lines lt 1 lw 2 lc rgb "#1080FF" title "VMAT/C-ion rectum", \
     "RBEmin_scan_cion_rectum.dat" using 1:3 with lines lt 3 lw 1 lc rgb "#1080FF" notitle, \
     "RBEmin_scan_cion_rectum.dat" using 1:5 with lines lt 3 lw 1 lc rgb "#1080FF" notitle
set key at 0.75,3
set xlabel "{/Symbol a}  [Gy^{-1}]"
set xrange [0:0.8]
set yrange [0:3.2]
plot "alpha_scan_proton_bladder.dat" using 1:2 with lines lt 1 lw 2 lc rgb "#CE0000" title "VMAT/IMPT bladder", \
     "alpha_scan_proton_bladder.dat" using 1:3 with lines lt 3 lw 1 lc rgb "#CE0000" notitle, \
     "alpha_scan_proton_bladder.dat" using 1:5 with lines lt 3 lw 1 lc rgb "#CE0000" notitle, \
     "alpha_scan_proton_rectum.dat" using 1:2 with lines lt 1 lw 2 lc rgb "#1080FF" title "VMAT/IMPT rectum", \
     "alpha_scan_proton_rectum.dat" using 1:3 with lines lt 3 lw 1 lc rgb "#1080FF" notitle, \
     "alpha_scan_proton_rectum.dat" using 1:5 with lines lt 3 lw 1 lc rgb "#1080FF" notitle
unset key
set xlabel "{/Symbol b}  [Gy^{-2}]"
set xrange [0.01:0.05]
set yrange [0:3.2]
plot "beta_scan_proton_bladder.dat" using 1:2 with lines lt 1 lw 2 lc rgb "#CE0000" title "VMAT/IMPT bladder", \
     "beta_scan_proton_bladder.dat" using 1:3 with lines lt 3 lw 1 lc rgb "#CE0000" notitle, \
     "beta_scan_proton_bladder.dat" using 1:5 with lines lt 3 lw 1 lc rgb "#CE0000" notitle, \
     "beta_scan_proton_rectum.dat" using 1:2 with lines lt 1 lw 2 lc rgb "#1080FF" title "VMAT/IMPT rectum", \
     "beta_scan_proton_rectum.dat" using 1:3 with lines lt 3 lw 1 lc rgb "#1080FF" notitle, \
     "beta_scan_proton_rectum.dat" using 1:5 with lines lt 3 lw 1 lc rgb "#1080FF" notitle
set xlabel "RBE_{max}"
set xrange [1.05:1.6]
set yrange [0:3.2]
plot "RBEmax_scan_proton_bladder.dat" using 1:2 with lines lt 1 lw 2 lc rgb "#CE0000" title "VMAT/IMPT bladder", \
     "RBEmax_scan_proton_bladder.dat" using 1:3 with lines lt 3 lw 1 lc rgb "#CE0000" notitle, \
     "RBEmax_scan_proton_bladder.dat" using 1:5 with lines lt 3 lw 1 lc rgb "#CE0000" notitle, \
     "RBEmax_scan_proton_rectum.dat" using 1:2 with lines lt 1 lw 2 lc rgb "#1080FF" title "VMAT/IMPT rectum", \
     "RBEmax_scan_proton_rectum.dat" using 1:3 with lines lt 3 lw 1 lc rgb "#1080FF" notitle, \
     "RBEmax_scan_proton_rectum.dat" using 1:5 with lines lt 3 lw 1 lc rgb "#1080FF" notitle
set xlabel "RBE_{min}"
set xrange [1:1.1]
set yrange [0:3.2]
plot "RBEmin_scan_proton_bladder.dat" using 1:2 with lines lt 1 lw 2 lc rgb "#CE0000" title "VMAT/IMPT bladder", \
     "RBEmin_scan_proton_bladder.dat" using 1:3 with lines lt 3 lw 1 lc rgb "#CE0000" notitle, \
     "RBEmin_scan_proton_bladder.dat" using 1:5 with lines lt 3 lw 1 lc rgb "#CE0000" notitle, \
     "RBEmin_scan_proton_rectum.dat" using 1:2 with lines lt 1 lw 2 lc rgb "#1080FF" title "VMAT/IMPT rectum", \
     "RBEmin_scan_proton_rectum.dat" using 1:3 with lines lt 3 lw 1 lc rgb "#1080FF" notitle, \
     "RBEmin_scan_proton_rectum.dat" using 1:5 with lines lt 3 lw 1 lc rgb "#1080FF" notitle
endef

export gnuplot_1d_summary_plot_script

summary_1d_plots.gp: scan_1D_plots.makefile
	echo "$$gnuplot_1d_summary_plot_script" > $@

summary_1d_plots.eps: summary_1d_plots.gp $(PLOT_FILES)
	gnuplot $<


define gnuplot_cion_1d_summary_plot_script
set terminal postscript eps color enhanced dashed size 16cm,12cm
set output "summary_cion_1d_plots.eps"
set multiplot layout 2,2 columnsfirst downwards
set ylabel "Relative Risk"
set xlabel "{/Symbol a}  [Gy^{-1}]"
set xrange [0:0.8]
set yrange [0:3.2]
set key at 0.75,3
plot "alpha_scan_cion_bladder.dat" using 1:2 with lines lt 1 lw 2 lc rgb "#CE0000" title "VMAT/C-ion bladder", \
     "alpha_scan_cion_bladder.dat" using 1:3 with lines lt 3 lw 1 lc rgb "#CE0000" notitle, \
     "alpha_scan_cion_bladder.dat" using 1:5 with lines lt 3 lw 1 lc rgb "#CE0000" notitle, \
     "alpha_scan_cion_rectum.dat" using 1:2 with lines lt 1 lw 2 lc rgb "#1080FF" title "VMAT/C-ion rectum", \
     "alpha_scan_cion_rectum.dat" using 1:3 with lines lt 3 lw 1 lc rgb "#1080FF" notitle, \
     "alpha_scan_cion_rectum.dat" using 1:5 with lines lt 3 lw 1 lc rgb "#1080FF" notitle
unset key
set xlabel "RBE_{max}"
set xrange [1:9]
set yrange [0:3.2]
plot "RBEmax_scan_cion_bladder.dat" using 1:2 with lines lt 1 lw 2 lc rgb "#CE0000" title "VMAT/C-ion bladder", \
     "RBEmax_scan_cion_bladder.dat" using 1:3 with lines lt 3 lw 1 lc rgb "#CE0000" notitle, \
     "RBEmax_scan_cion_bladder.dat" using 1:5 with lines lt 3 lw 1 lc rgb "#CE0000" notitle, \
     "RBEmax_scan_cion_rectum.dat" using 1:2 with lines lt 1 lw 2 lc rgb "#1080FF" title "VMAT/C-ion rectum", \
     "RBEmax_scan_cion_rectum.dat" using 1:3 with lines lt 3 lw 1 lc rgb "#1080FF" notitle, \
     "RBEmax_scan_cion_rectum.dat" using 1:5 with lines lt 3 lw 1 lc rgb "#1080FF" notitle
set xlabel "{/Symbol b}  [Gy^{-2}]"
set xrange [0.01:0.05]
set yrange [0:3.2]
plot "beta_scan_cion_bladder.dat" using 1:2 with lines lt 1 lw 2 lc rgb "#CE0000" title "VMAT/C-ion bladder", \
     "beta_scan_cion_bladder.dat" using 1:3 with lines lt 3 lw 1 lc rgb "#CE0000" notitle, \
     "beta_scan_cion_bladder.dat" using 1:5 with lines lt 3 lw 1 lc rgb "#CE0000" notitle, \
     "beta_scan_cion_rectum.dat" using 1:2 with lines lt 1 lw 2 lc rgb "#1080FF" title "VMAT/C-ion rectum", \
     "beta_scan_cion_rectum.dat" using 1:3 with lines lt 3 lw 1 lc rgb "#1080FF" notitle, \
     "beta_scan_cion_rectum.dat" using 1:5 with lines lt 3 lw 1 lc rgb "#1080FF" notitle
set xlabel "RBE_{min}"
set xrange [1:1.4]
set yrange [0:3.2]
plot "RBEmin_scan_cion_bladder.dat" using 1:2 with lines lt 1 lw 2 lc rgb "#CE0000" title "VMAT/C-ion bladder", \
     "RBEmin_scan_cion_bladder.dat" using 1:3 with lines lt 3 lw 1 lc rgb "#CE0000" notitle, \
     "RBEmin_scan_cion_bladder.dat" using 1:5 with lines lt 3 lw 1 lc rgb "#CE0000" notitle, \
     "RBEmin_scan_cion_rectum.dat" using 1:2 with lines lt 1 lw 2 lc rgb "#1080FF" title "VMAT/C-ion rectum", \
     "RBEmin_scan_cion_rectum.dat" using 1:3 with lines lt 3 lw 1 lc rgb "#1080FF" notitle, \
     "RBEmin_scan_cion_rectum.dat" using 1:5 with lines lt 3 lw 1 lc rgb "#1080FF" notitle
endef


define gnuplot_proton_1d_summary_plot_script
set terminal postscript eps color enhanced dashed size 16cm,12cm
set output "summary_proton_1d_plots.eps"
set multiplot layout 2,2 columnsfirst downwards
set ylabel "Relative Risk"
set xlabel "{/Symbol a}  [Gy^{-1}]"
set xrange [0:0.8]
set yrange [0:3.2]
set key at 0.75,3
plot "alpha_scan_proton_bladder.dat" using 1:2 with lines lt 1 lw 2 lc rgb "#CE0000" title "VMAT/IMPT bladder", \
     "alpha_scan_proton_bladder.dat" using 1:3 with lines lt 3 lw 1 lc rgb "#CE0000" notitle, \
     "alpha_scan_proton_bladder.dat" using 1:5 with lines lt 3 lw 1 lc rgb "#CE0000" notitle, \
     "alpha_scan_proton_rectum.dat" using 1:2 with lines lt 1 lw 2 lc rgb "#1080FF" title "VMAT/IMPT rectum", \
     "alpha_scan_proton_rectum.dat" using 1:3 with lines lt 3 lw 1 lc rgb "#1080FF" notitle, \
     "alpha_scan_proton_rectum.dat" using 1:5 with lines lt 3 lw 1 lc rgb "#1080FF" notitle
unset key
set xlabel "RBE_{max}"
set xrange [1.05:1.6]
set yrange [0:3.2]
plot "RBEmax_scan_proton_bladder.dat" using 1:2 with lines lt 1 lw 2 lc rgb "#CE0000" title "VMAT/IMPT bladder", \
     "RBEmax_scan_proton_bladder.dat" using 1:3 with lines lt 3 lw 1 lc rgb "#CE0000" notitle, \
     "RBEmax_scan_proton_bladder.dat" using 1:5 with lines lt 3 lw 1 lc rgb "#CE0000" notitle, \
     "RBEmax_scan_proton_rectum.dat" using 1:2 with lines lt 1 lw 2 lc rgb "#1080FF" title "VMAT/IMPT rectum", \
     "RBEmax_scan_proton_rectum.dat" using 1:3 with lines lt 3 lw 1 lc rgb "#1080FF" notitle, \
     "RBEmax_scan_proton_rectum.dat" using 1:5 with lines lt 3 lw 1 lc rgb "#1080FF" notitle
set xlabel "{/Symbol b}  [Gy^{-2}]"
set xrange [0.01:0.05]
set yrange [0:3.2]
plot "beta_scan_proton_bladder.dat" using 1:2 with lines lt 1 lw 2 lc rgb "#CE0000" title "VMAT/IMPT bladder", \
     "beta_scan_proton_bladder.dat" using 1:3 with lines lt 3 lw 1 lc rgb "#CE0000" notitle, \
     "beta_scan_proton_bladder.dat" using 1:5 with lines lt 3 lw 1 lc rgb "#CE0000" notitle, \
     "beta_scan_proton_rectum.dat" using 1:2 with lines lt 1 lw 2 lc rgb "#1080FF" title "VMAT/IMPT rectum", \
     "beta_scan_proton_rectum.dat" using 1:3 with lines lt 3 lw 1 lc rgb "#1080FF" notitle, \
     "beta_scan_proton_rectum.dat" using 1:5 with lines lt 3 lw 1 lc rgb "#1080FF" notitle
set xlabel "RBE_{min}"
set xrange [1:1.1]
set yrange [0:3.2]
plot "RBEmin_scan_proton_bladder.dat" using 1:2 with lines lt 1 lw 2 lc rgb "#CE0000" title "VMAT/IMPT bladder", \
     "RBEmin_scan_proton_bladder.dat" using 1:3 with lines lt 3 lw 1 lc rgb "#CE0000" notitle, \
     "RBEmin_scan_proton_bladder.dat" using 1:5 with lines lt 3 lw 1 lc rgb "#CE0000" notitle, \
     "RBEmin_scan_proton_rectum.dat" using 1:2 with lines lt 1 lw 2 lc rgb "#1080FF" title "VMAT/IMPT rectum", \
     "RBEmin_scan_proton_rectum.dat" using 1:3 with lines lt 3 lw 1 lc rgb "#1080FF" notitle, \
     "RBEmin_scan_proton_rectum.dat" using 1:5 with lines lt 3 lw 1 lc rgb "#1080FF" notitle
endef


export gnuplot_cion_1d_summary_plot_script
export gnuplot_proton_1d_summary_plot_script

summary_cion_1d_plots.gp: scan_1D_plots.makefile
	echo "$$gnuplot_cion_1d_summary_plot_script" > $@

summary_proton_1d_plots.gp: scan_1D_plots.makefile
	echo "$$gnuplot_proton_1d_summary_plot_script" > $@

summary_cion_1d_plots.eps: summary_cion_1d_plots.gp $(PLOT_FILES)
	gnuplot $<

summary_proton_1d_plots.eps: summary_proton_1d_plots.gp $(PLOT_FILES)
	gnuplot $<
