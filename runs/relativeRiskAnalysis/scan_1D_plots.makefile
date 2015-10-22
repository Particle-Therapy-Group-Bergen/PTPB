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

all: $(PLOT_FILES)

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
	./make1DscanPlots.sh $< $(subst .eps,.dat,$@) '1:0.1:9' 5 $@ 'RBE_{min}  [Gy]'

RBEmax_scan_proton_%.eps: RBEmax_scan_proton_%_sample_data.mat scan_1D_plots.makefile make1DscanPlots.sh
	./make1DscanPlots.sh $< $(subst .eps,.dat,$@) '1.05:0.01:1.6' 5 $@ 'RBE_{min}  [Gy]'

