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

# This generates the scan plots after hypercubes.makefile has finished creating
# the hypercube files.

FILEIDS = $(foreach N,1 2 3 4,$(foreach M,1 2 3 4,$(N)$(M)))

FILETYPES = contour_labels_FILEID.dat \
            contours_FILEID.dat \
            contours_high_FILEID.dat \
            contours_low_FILEID.dat \
            processed_contours_FILEID.dat \
            relative_risk_FILEID.dat \
            relative_risk_high_FILEID.dat \
            relative_risk_low_FILEID.dat

PLOT_FILES = $(foreach N,cion_bladder cion_rectum proton_bladder proton_rectum,$(N).eps)

OUTPUT_FILES = $(foreach N,$(FILEIDS),$(subst FILEID,cion_bladder_$(N),$(FILETYPES))) \
               $(foreach N,$(FILEIDS),$(subst FILEID,cion_rectum_$(N),$(FILETYPES))) \
               $(foreach N,$(FILEIDS),$(subst FILEID,proton_bladder_$(N),$(FILETYPES))) \
               $(foreach N,$(FILEIDS),$(subst FILEID,proton_rectum_$(N),$(FILETYPES)))

.PHONY: all clean cleanall

all: $(PLOT_FILES)

clean:
	rm -rf $(OUTPUT_FILES)

cleanall: clean
	rm -rf $(PLOT_FILES)

%.eps: hypercube_%.mat makeScanPlot.sh makeContourLabels.py
	./makeScanPlot.sh $*
