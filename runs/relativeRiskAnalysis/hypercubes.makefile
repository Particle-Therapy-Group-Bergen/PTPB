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

# This generates the hypercube files with scans of Relative Risk over the alpha,
# beta and RBE min/max parameters.

OUTPUT_FILES = hypercube_cion_bladder.mat \
               hypercube_cion_rectum.mat \
               hypercube_proton_bladder.mat \
               hypercube_proton_rectum.mat

OCTAVE = octave -q --path ~/bin/PTPB_mfiles

.PHONY: all clean cleanall

all: $(OUTPUT_FILES)

clean:
	rm -rf $(foreach N,$(OUTPUT_FILES),log_$(N).txt)

cleanall: clean
	rm -rf make_hypercube $(OUTPUT_FILES)


hypercube_cion_bladder.mat: make_hypercube
	./make_hypercube $@ 'data/VMATdvh/vmat%d.mat' 25 1 \
		'data/CionDataPhysicalDose/HUH%dphysical_dvh.mat' 12 1 \
		Bladder '[1.05 1.15 1.25 1.35]' '[2 4 6 8]' > log_$@.txt 2>&1

hypercube_cion_rectum.mat: make_hypercube
	./make_hypercube $@ 'data/VMATdvh/vmat%d.mat' 25 1 \
		'data/CionDataPhysicalDose/HUH%dphysical_dvh.mat' 12 1 \
		Rectum '[1.05 1.15 1.25 1.35]' '[2 4 6 8]' > log_$@.txt 2>&1

hypercube_proton_bladder.mat: make_hypercube
	./make_hypercube $@ 'data/VMATdvh/vmat%d.mat' 25 1 \
		'data/IMPTdvh/impt%d.mat' 25 1/1.1 \
		Bladder '[1.00 1.03 1.06 1.09]' '[1.10 1.25 1.40 1.55]' > log_$@.txt 2>&1

hypercube_proton_rectum.mat: make_hypercube
	./make_hypercube $@ 'data/VMATdvh/vmat%d.mat' 25 1 \
		'data/IMPTdvh/impt%d.mat' 25 1/1.1 \
		Rectum '[1.00 1.03 1.06 1.09]' '[1.10 1.25 1.40 1.55]' > log_$@.txt 2>&1


define make_hypercube
#!/bin/sh
exec $(OCTAVE) <<EOF
patients = [12, 33, 35, 36, 37, 39, 41, 42, 43, 44];
filepat1 = '$$2';
n1 = $$3;
scale1 = $$4;
filepat2 = '$$5';
n2 = $$6;
scale2 = $$7;
organ = '$$8';
RBEmin_vec = $$9;
RBEmax_vec = $$10;
[alpha, beta, RBEmin, RBEmax, RR, RRerr] = makeRelativeRiskHyperCube(
                   organ, patients, filepat1, filepat2, n1, n2, scale1, scale2,
                   RBEmin_vec, RBEmax_vec);
save('-v7', '$$1', 'alpha', 'beta', 'RBEmin', 'RBEmax', 'RR', 'RRerr');
EOF
endef
export make_hypercube

make_hypercube: hypercubes.makefile
	echo "$$make_hypercube" > $@
	chmod +x $@
