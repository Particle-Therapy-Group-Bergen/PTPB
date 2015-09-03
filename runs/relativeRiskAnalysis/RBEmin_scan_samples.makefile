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

# This makefile is used for producing 1D scans of the RBE_max parameter
# around the nominal point.

ITERATIONS = 10
N_SAMPLES = 1000

ALPHA_CION_BLADDER_MEAN = 0.25
ALPHA_CION_BLADDER_STD = 0.075
BETA_CION_BLADDER_MEAN = 0.033
BETA_CION_BLADDER_STD = 0.0055
RBEMIN_CION_BLADDER_LOW = 1.0
RBEMIN_CION_BLADDER_HIGH = 1.4
RBEMAX_CION_BLADDER_LOW = 5
RBEMAX_CION_BLADDER_MODE = 6
RBEMAX_CION_BLADDER_HIGH = 7

ALPHA_CION_RECTUM_MEAN = $(ALPHA_CION_BLADDER_MEAN)
ALPHA_CION_RECTUM_STD = $(ALPHA_CION_BLADDER_STD)
BETA_CION_RECTUM_MEAN = 0.046
BETA_CION_RECTUM_STD = 0.0077
RBEMIN_CION_RECTUM_LOW = $(RBEMIN_CION_BLADDER_LOW)
RBEMIN_CION_RECTUM_HIGH = $(RBEMIN_CION_BLADDER_HIGH)
RBEMAX_CION_RECTUM_LOW = $(RBEMAX_CION_BLADDER_LOW)
RBEMAX_CION_RECTUM_MODE = $(RBEMAX_CION_BLADDER_MODE)
RBEMAX_CION_RECTUM_HIGH = $(RBEMAX_CION_BLADDER_HIGH)

ALPHA_PROTON_BLADDER_MEAN = $(ALPHA_CION_BLADDER_MEAN)
ALPHA_PROTON_BLADDER_STD = $(ALPHA_CION_BLADDER_STD)
BETA_PROTON_BLADDER_MEAN = $(BETA_CION_BLADDER_MEAN)
BETA_PROTON_BLADDER_STD = $(BETA_CION_BLADDER_STD)
RBEMIN_PROTON_BLADDER_LOW = 1.00
RBEMIN_PROTON_BLADDER_HIGH = 1.1
RBEMAX_PROTON_BLADDER_LOW = 1.2
RBEMAX_PROTON_BLADDER_MODE = 1.25
RBEMAX_PROTON_BLADDER_HIGH = 1.3

ALPHA_PROTON_RECTUM_MEAN = $(ALPHA_CION_RECTUM_MEAN)
ALPHA_PROTON_RECTUM_STD = $(ALPHA_CION_RECTUM_STD)
BETA_PROTON_RECTUM_MEAN = $(BETA_CION_RECTUM_MEAN)
BETA_PROTON_RECTUM_STD = $(BETA_CION_RECTUM_STD)
RBEMIN_PROTON_RECTUM_LOW = $(RBEMIN_PROTON_BLADDER_LOW)
RBEMIN_PROTON_RECTUM_HIGH = $(RBEMIN_PROTON_BLADDER_HIGH)
RBEMAX_PROTON_RECTUM_LOW = $(RBEMAX_PROTON_BLADDER_LOW)
RBEMAX_PROTON_RECTUM_MODE = $(RBEMAX_PROTON_BLADDER_MODE)
RBEMAX_PROTON_RECTUM_HIGH = $(RBEMAX_PROTON_BLADDER_HIGH)

OCTAVE = octave -q --path ~/bin/PTPB_mfiles

FILEIDS = $(shell N=1; while test $$N -le $(ITERATIONS); do echo $$N ; N=$$((N+1)) ; done)

CION_BLADDER_FILELIST = $(foreach N,$(FILEIDS),RBEmin_scan_cion_bladder_samples_$(N).mat)
CION_RECTUM_FILELIST = $(foreach N,$(FILEIDS),RBEmin_scan_cion_rectum_samples_$(N).mat)
PROTON_BLADDER_FILELIST = $(foreach N,$(FILEIDS),RBEmin_scan_proton_bladder_samples_$(N).mat)
PROTON_RECTUM_FILELIST = $(foreach N,$(FILEIDS),RBEmin_scan_proton_rectum_samples_$(N).mat)
FILELIST = $(CION_BLADDER_FILELIST) $(CION_RECTUM_FILELIST) $(PROTON_BLADDER_FILELIST) $(PROTON_RECTUM_FILELIST)

OUTPUT_FILES = RBEmin_scan_cion_bladder_sample_data.mat \
               RBEmin_scan_cion_rectum_sample_data.mat \
               RBEmin_scan_proton_bladder_sample_data.mat \
               RBEmin_scan_proton_rectum_sample_data.mat


.PHONY: all clean

all: $(OUTPUT_FILES)

clean:
	rm -rf $(FILELIST) $(foreach N,$(FILELIST),log_$(N).txt)

cleanall: clean
	rm -rf make_RBEmin_scan_samples merge_RBEmin_scan_samples $(OUTPUT_FILES)


define merge_RBEmin_scan_samples
#!/bin/sh
outfile="$$1"
shift
if test -f "$$outfile" ; then
    mv "$$outfile" "$${outfile}.backup"
    files="'$${outfile}.backup', "
else
    files=""
fi
for N ; do
    files="$$files'$$N', "
done
exec $(OCTAVE) <<EOF
files = {$$files};
Results = [];
for n = 1:length(files)
    data = load(files{n}, 'Results').Results;
    Results = [Results ; data];
end
save('-v7', '$$outfile', 'Results');
EOF
endef
export merge_RBEmin_scan_samples

merge_RBEmin_scan_samples: RBEmin_scan_samples.makefile
	echo "$$merge_RBEmin_scan_samples" > $@
	chmod +x $@


define make_RBEmin_scan_samples
#!/bin/sh
exec $(OCTAVE) <<EOF
Nsamples = $$2;
patients = [12, 33, 35, 36, 37, 39, 41, 42, 43, 44];
filepat1 = '$$3';
n1 = $$4;
scale1 = $$5;
filepat2 = '$$6';
n2 = $$7;
scale2 = $$8;
organ = '$$9';
alpha_distrib = struct('type', 'gaus', 'params', {{$$10, $$11}});
beta_distrib = struct('type', 'gaus', 'params', {{$$12, $$13}});
RBEmin_distrib = struct('type', 'box', 'params', {{$$14, $$15}});
RBEmax_distrib = struct('type', 'triangle', 'params', {{$$16, $$17, $$18}});
opts = struct('integration_method', 'trapz',
              'integration_tolerance', 1e-4,
              'interpolation_method', 'pchip',
              'sample_dvh', 1,
              'bootstrap_samples', 10,
              'bootstrap_method', 'random');
namemap = {
        'Bladder_P', 'Bladder';
        'Rectum_P_MT', 'Rectum';
    };
Results = sampleMeanRelativeRisk(Nsamples, filepat1, filepat2, patients, organ,
                           n1, n2, scale1, scale2, alpha_distrib, beta_distrib,
                           RBEmin_distrib, RBEmax_distrib, opts, namemap);
save('-v7', '$$1', 'Results');
EOF
endef
export make_RBEmin_scan_samples

make_RBEmin_scan_samples: RBEmin_scan_samples.makefile
	echo "$$make_RBEmin_scan_samples" > $@
	chmod +x $@


$(CION_BLADDER_FILELIST): make_RBEmin_scan_samples
	./make_RBEmin_scan_samples $@ $(N_SAMPLES) 'data/VMATdvh/vmat%d.mat' 25 1 \
		'data/CionDataPhysicalDose/HUH%dphysical_dvh.mat' 12 1 Bladder \
		$(ALPHA_CION_BLADDER_MEAN) $(ALPHA_CION_BLADDER_STD) \
		$(BETA_CION_BLADDER_MEAN) $(BETA_CION_BLADDER_STD) \
		$(RBEMIN_CION_BLADDER_LOW) $(RBEMIN_CION_BLADDER_HIGH) \
		$(RBEMAX_CION_BLADDER_LOW) $(RBEMAX_CION_BLADDER_MODE) $(RBEMAX_CION_BLADDER_HIGH) \
		> log_$@.txt 2>&1


RBEmin_scan_cion_bladder_sample_data.mat: merge_RBEmin_scan_samples $(CION_BLADDER_FILELIST)
	./merge_RBEmin_scan_samples $@ $(CION_BLADDER_FILELIST) && \
		rm -f $@.backup $(CION_BLADDER_FILELIST) $(foreach N,$(CION_BLADDER_FILELIST),log_$(N).txt)


$(CION_RECTUM_FILELIST): make_RBEmin_scan_samples
	./make_RBEmin_scan_samples $@ $(N_SAMPLES) 'data/VMATdvh/vmat%d.mat' 25 1 \
		'data/CionDataPhysicalDose/HUH%dphysical_dvh.mat' 12 1 Rectum \
		$(ALPHA_CION_RECTUM_MEAN) $(ALPHA_CION_RECTUM_STD) \
		$(BETA_CION_RECTUM_MEAN) $(BETA_CION_RECTUM_STD) \
		$(RBEMIN_CION_RECTUM_LOW) $(RBEMIN_CION_RECTUM_HIGH) \
		$(RBEMAX_CION_RECTUM_LOW) $(RBEMAX_CION_RECTUM_MODE) $(RBEMAX_CION_RECTUM_HIGH) \
		> log_$@.txt 2>&1


RBEmin_scan_cion_rectum_sample_data.mat: merge_RBEmin_scan_samples $(CION_RECTUM_FILELIST)
	./merge_RBEmin_scan_samples $@ $(CION_RECTUM_FILELIST) && \
		rm -f $@.backup $(CION_RECTUM_FILELIST) $(foreach N,$(CION_RECTUM_FILELIST),log_$(N).txt)


$(PROTON_BLADDER_FILELIST): make_RBEmin_scan_samples
	./make_RBEmin_scan_samples $@ $(N_SAMPLES) 'data/VMATdvh/vmat%d.mat' 25 1 \
		'data/IMPTdvh/impt%d.mat' 25 1 Bladder \
		$(ALPHA_PROTON_BLADDER_MEAN) $(ALPHA_PROTON_BLADDER_STD) \
		$(BETA_PROTON_BLADDER_MEAN) $(BETA_PROTON_BLADDER_STD) \
		$(RBEMIN_PROTON_BLADDER_LOW) $(RBEMIN_PROTON_BLADDER_HIGH) \
		$(RBEMAX_PROTON_BLADDER_LOW) $(RBEMAX_PROTON_BLADDER_MODE) $(RBEMAX_PROTON_BLADDER_HIGH) \
		> log_$@.txt 2>&1


RBEmin_scan_proton_bladder_sample_data.mat: merge_RBEmin_scan_samples $(PROTON_BLADDER_FILELIST)
	./merge_RBEmin_scan_samples $@ $(PROTON_BLADDER_FILELIST) && \
		rm -f $@.backup $(PROTON_BLADDER_FILELIST) $(foreach N,$(PROTON_BLADDER_FILELIST),log_$(N).txt)


$(PROTON_RECTUM_FILELIST): make_RBEmin_scan_samples
	./make_RBEmin_scan_samples $@ $(N_SAMPLES) 'data/VMATdvh/vmat%d.mat' 25 1 \
		'data/IMPTdvh/impt%d.mat' 25 1 Rectum \
		$(ALPHA_PROTON_RECTUM_MEAN) $(ALPHA_PROTON_RECTUM_STD) \
		$(BETA_PROTON_RECTUM_MEAN) $(BETA_PROTON_RECTUM_STD) \
		$(RBEMIN_PROTON_RECTUM_LOW) $(RBEMIN_PROTON_RECTUM_HIGH) \
		$(RBEMAX_PROTON_RECTUM_LOW) $(RBEMAX_PROTON_RECTUM_MODE) $(RBEMAX_PROTON_RECTUM_HIGH) \
		> log_$@.txt 2>&1


RBEmin_scan_proton_rectum_sample_data.mat: merge_RBEmin_scan_samples $(PROTON_RECTUM_FILELIST)
	./merge_RBEmin_scan_samples $@ $(PROTON_RECTUM_FILELIST) && \
		rm -f $@.backup $(PROTON_RECTUM_FILELIST) $(foreach N,$(PROTON_RECTUM_FILELIST),log_$(N).txt)

