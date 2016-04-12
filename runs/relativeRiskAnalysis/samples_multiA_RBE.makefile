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

# This makefile will generate samples for nominal C-ion and Proton Relative
# Risk values.

ITERATIONS = 100
N_SAMPLES = 1000

ALPHA_CION_BLADDER_MEAN = 0.25
ALPHA_CION_BLADDER_STD = 0.075
BETA_CION_BLADDER_MEAN = 0.033
BETA_CION_BLADDER_STD = 0.0055
RBEMIN1_CION_BLADDER_LOW = 1.2
RBEMIN1_CION_BLADDER_MODE = 1.25
RBEMIN1_CION_BLADDER_HIGH = 1.3
RBEMAX1_CION_BLADDER_LOW = 7.5
RBEMAX1_CION_BLADDER_MODE = 9
RBEMAX1_CION_BLADDER_HIGH = 10.5
RBEMIN2_CION_BLADDER_LOW = $(RBEMIN1_CION_BLADDER_LOW)
RBEMIN2_CION_BLADDER_MODE = $(RBEMIN1_CION_BLADDER_MODE)
RBEMIN2_CION_BLADDER_HIGH = $(RBEMIN1_CION_BLADDER_HIGH)
RBEMAX2_CION_BLADDER_LOW = 5
RBEMAX2_CION_BLADDER_MODE = 6
RBEMAX2_CION_BLADDER_HIGH = 7

ALPHA_CION_RECTUM_MEAN = $(ALPHA_CION_BLADDER_MEAN)
ALPHA_CION_RECTUM_STD = $(ALPHA_CION_BLADDER_STD)
BETA_CION_RECTUM_MEAN = 0.046
BETA_CION_RECTUM_STD = 0.0077
RBEMIN1_CION_RECTUM_LOW = $(RBEMIN1_CION_BLADDER_LOW)
RBEMIN1_CION_RECTUM_MODE = $(RBEMIN1_CION_BLADDER_MODE)
RBEMIN1_CION_RECTUM_HIGH = $(RBEMIN1_CION_BLADDER_HIGH)
RBEMAX1_CION_RECTUM_LOW = $(RBEMAX1_CION_BLADDER_LOW)
RBEMAX1_CION_RECTUM_MODE = $(RBEMAX1_CION_BLADDER_MODE)
RBEMAX1_CION_RECTUM_HIGH = $(RBEMAX1_CION_BLADDER_HIGH)
RBEMIN2_CION_RECTUM_LOW = $(RBEMIN2_CION_BLADDER_LOW)
RBEMIN2_CION_RECTUM_MODE = $(RBEMIN2_CION_BLADDER_MODE)
RBEMIN2_CION_RECTUM_HIGH = $(RBEMIN2_CION_BLADDER_HIGH)
RBEMAX2_CION_RECTUM_LOW = $(RBEMAX2_CION_BLADDER_LOW)
RBEMAX2_CION_RECTUM_MODE = $(RBEMAX2_CION_BLADDER_MODE)
RBEMAX2_CION_RECTUM_HIGH = $(RBEMAX2_CION_BLADDER_HIGH)

ALPHA_PROTON_BLADDER_MEAN = $(ALPHA_CION_BLADDER_MEAN)
ALPHA_PROTON_BLADDER_STD = $(ALPHA_CION_BLADDER_STD)
BETA_PROTON_BLADDER_MEAN = $(BETA_CION_BLADDER_MEAN)
BETA_PROTON_BLADDER_STD = $(BETA_CION_BLADDER_STD)
RBEMIN1_PROTON_BLADDER_LOW = 1.01
RBEMIN1_PROTON_BLADDER_MODE = 1.03
RBEMIN1_PROTON_BLADDER_HIGH = 1.05
RBEMAX1_PROTON_BLADDER_LOW = 1.8
RBEMAX1_PROTON_BLADDER_MODE = 1.875
RBEMAX1_PROTON_BLADDER_HIGH = 1.95
RBEMIN2_PROTON_BLADDER_LOW = $(RBEMIN1_PROTON_BLADDER_LOW)
RBEMIN2_PROTON_BLADDER_MODE = $(RBEMIN1_PROTON_BLADDER_MODE)
RBEMIN2_PROTON_BLADDER_HIGH = $(RBEMIN1_PROTON_BLADDER_HIGH)
RBEMAX2_PROTON_BLADDER_LOW = 1.2
RBEMAX2_PROTON_BLADDER_MODE = 1.25
RBEMAX2_PROTON_BLADDER_HIGH = 1.3

ALPHA_PROTON_RECTUM_MEAN = $(ALPHA_CION_RECTUM_MEAN)
ALPHA_PROTON_RECTUM_STD = $(ALPHA_CION_RECTUM_STD)
BETA_PROTON_RECTUM_MEAN = $(BETA_CION_RECTUM_MEAN)
BETA_PROTON_RECTUM_STD = $(BETA_CION_RECTUM_STD)
RBEMIN1_PROTON_RECTUM_LOW = $(RBEMIN1_PROTON_BLADDER_LOW)
RBEMIN1_PROTON_RECTUM_MODE = $(RBEMIN1_PROTON_BLADDER_MODE)
RBEMIN1_PROTON_RECTUM_HIGH = $(RBEMIN1_PROTON_BLADDER_HIGH)
RBEMAX1_PROTON_RECTUM_LOW = $(RBEMAX1_PROTON_BLADDER_LOW)
RBEMAX1_PROTON_RECTUM_MODE = $(RBEMAX1_PROTON_BLADDER_MODE)
RBEMAX1_PROTON_RECTUM_HIGH = $(RBEMAX1_PROTON_BLADDER_HIGH)
RBEMIN2_PROTON_RECTUM_LOW = $(RBEMIN2_PROTON_BLADDER_LOW)
RBEMIN2_PROTON_RECTUM_MODE = $(RBEMIN2_PROTON_BLADDER_MODE)
RBEMIN2_PROTON_RECTUM_HIGH = $(RBEMIN2_PROTON_BLADDER_HIGH)
RBEMAX2_PROTON_RECTUM_LOW = $(RBEMAX2_PROTON_BLADDER_LOW)
RBEMAX2_PROTON_RECTUM_MODE = $(RBEMAX2_PROTON_BLADDER_MODE)
RBEMAX2_PROTON_RECTUM_HIGH = $(RBEMAX2_PROTON_BLADDER_HIGH)

OCTAVE = octave -q --path ~/bin/PTPB_mfiles

FILEIDS = $(shell N=1; while test $$N -le $(ITERATIONS); do echo $$N ; N=$$((N+1)) ; done)

CION_BLADDER_FILELIST = $(foreach N,$(FILEIDS),cion_bladder_samples_multiA_RBE_$(N).mat)
CION_RECTUM_FILELIST = $(foreach N,$(FILEIDS),cion_rectum_samples_multiA_RBE_$(N).mat)
PROTON_BLADDER_FILELIST = $(foreach N,$(FILEIDS),proton_bladder_samples_multiA_RBE_$(N).mat)
PROTON_RECTUM_FILELIST = $(foreach N,$(FILEIDS),proton_rectum_samples_multiA_RBE_$(N).mat)
FILELIST = $(CION_BLADDER_FILELIST) $(CION_RECTUM_FILELIST) $(PROTON_BLADDER_FILELIST) $(PROTON_RECTUM_FILELIST)

OUTPUT_FILES = cion_bladder_sample_multiA_RBE_data.mat \
               cion_rectum_sample_multiA_RBE_data.mat \
               proton_bladder_sample_multiA_RBE_data.mat \
               proton_rectum_sample_multiA_RBE_data.mat

.PHONY: all clean

all: $(OUTPUT_FILES)

clean:
	rm -rf $(FILELIST) $(foreach N,$(FILELIST),log_$(N).txt)

cleanall: clean
	rm -rf make_samples_multiA_RBE merge_samples_multiA_RBE $(OUTPUT_FILES)


define merge_samples_multiA_RBE
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
export merge_samples_multiA_RBE

merge_samples_multiA_RBE: samples_multiA_RBE.makefile
	echo "$$merge_samples_multiA_RBE" > $@
	chmod +x $@


define make_samples_multiA_RBE
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
RBEmin1_distrib = struct('type', 'triangle', 'params', {{$$14, $$15, $$16}});
RBEmax1_distrib = struct('type', 'triangle', 'params', {{$$17, $$18, $$19}});
RBEmin2_distrib = struct('type', 'triangle', 'params', {{$$20, $$21, $$22}});
RBEmax2_distrib = struct('type', 'triangle', 'params', {{$$23, $$24, $$25}});
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
Results = sampleMeanRelativeRiskMultiRBE(Nsamples, filepat1, filepat2, patients,
                        organ, n1, n2, scale1, scale2,
                        alpha_distrib, beta_distrib,
                        RBEmin1_distrib, RBEmax1_distrib,
                        RBEmin2_distrib, RBEmax2_distrib,
                        opts, namemap);
save('-v7', '$$1', 'Results');
EOF
endef
export make_samples_multiA_RBE

make_samples_multiA_RBE: samples_multiA_RBE.makefile
	echo "$$make_samples_multiA_RBE" > $@
	chmod +x $@


$(CION_BLADDER_FILELIST): make_samples_multiA_RBE
	./make_samples_multiA_RBE $@ $(N_SAMPLES) 'data/VMATdvh/vmat%d.mat' 25 1 \
		'data/CionDataPhysicalDose/HUH%dphysical_dvh.mat' 12 1 Bladder \
		$(ALPHA_CION_BLADDER_MEAN) $(ALPHA_CION_BLADDER_STD) \
		$(BETA_CION_BLADDER_MEAN) $(BETA_CION_BLADDER_STD) \
		$(RBEMIN1_CION_BLADDER_LOW) $(RBEMIN1_CION_BLADDER_MODE) $(RBEMIN1_CION_BLADDER_HIGH) \
		$(RBEMAX1_CION_BLADDER_LOW) $(RBEMAX1_CION_BLADDER_MODE) $(RBEMAX1_CION_BLADDER_HIGH) \
		$(RBEMIN2_CION_BLADDER_LOW) $(RBEMIN2_CION_BLADDER_MODE) $(RBEMIN2_CION_BLADDER_HIGH) \
		$(RBEMAX2_CION_BLADDER_LOW) $(RBEMAX2_CION_BLADDER_MODE) $(RBEMAX2_CION_BLADDER_HIGH) \
		> log_$@.txt 2>&1


cion_bladder_sample_multiA_RBE_data.mat: merge_samples_multiA_RBE $(CION_BLADDER_FILELIST)
	./merge_samples_multiA_RBE $@ $(CION_BLADDER_FILELIST) && \
		rm -f $@.backup $(CION_BLADDER_FILELIST) $(foreach N,$(CION_BLADDER_FILELIST),log_$(N).txt)


$(CION_RECTUM_FILELIST): make_samples_multiA_RBE
	./make_samples_multiA_RBE $@ $(N_SAMPLES) 'data/VMATdvh/vmat%d.mat' 25 1 \
		'data/CionDataPhysicalDose/HUH%dphysical_dvh.mat' 12 1 Rectum \
		$(ALPHA_CION_RECTUM_MEAN) $(ALPHA_CION_RECTUM_STD) \
		$(BETA_CION_RECTUM_MEAN) $(BETA_CION_RECTUM_STD) \
		$(RBEMIN1_CION_RECTUM_LOW) $(RBEMIN1_CION_RECTUM_MODE) $(RBEMIN1_CION_RECTUM_HIGH) \
		$(RBEMAX1_CION_RECTUM_LOW) $(RBEMAX1_CION_RECTUM_MODE) $(RBEMAX1_CION_RECTUM_HIGH) \
		$(RBEMIN2_CION_RECTUM_LOW) $(RBEMIN2_CION_RECTUM_MODE) $(RBEMIN2_CION_RECTUM_HIGH) \
		$(RBEMAX2_CION_RECTUM_LOW) $(RBEMAX2_CION_RECTUM_MODE) $(RBEMAX2_CION_RECTUM_HIGH) \
		> log_$@.txt 2>&1


cion_rectum_sample_multiA_RBE_data.mat: merge_samples_multiA_RBE $(CION_RECTUM_FILELIST)
	./merge_samples_multiA_RBE $@ $(CION_RECTUM_FILELIST) && \
		rm -f $@.backup $(CION_RECTUM_FILELIST) $(foreach N,$(CION_RECTUM_FILELIST),log_$(N).txt)


$(PROTON_BLADDER_FILELIST): make_samples_multiA_RBE
	./make_samples_multiA_RBE $@ $(N_SAMPLES) 'data/VMATdvh/vmat%d.mat' 25 1 \
		'data/IMPTdvh/impt%d.mat' 25 1/1.1 Bladder \
		$(ALPHA_PROTON_BLADDER_MEAN) $(ALPHA_PROTON_BLADDER_STD) \
		$(BETA_PROTON_BLADDER_MEAN) $(BETA_PROTON_BLADDER_STD) \
		$(RBEMIN1_PROTON_BLADDER_LOW) $(RBEMIN1_PROTON_BLADDER_MODE) $(RBEMIN1_PROTON_BLADDER_HIGH) \
		$(RBEMAX1_PROTON_BLADDER_LOW) $(RBEMAX1_PROTON_BLADDER_MODE) $(RBEMAX1_PROTON_BLADDER_HIGH) \
		$(RBEMIN2_PROTON_BLADDER_LOW) $(RBEMIN2_PROTON_BLADDER_MODE) $(RBEMIN2_PROTON_BLADDER_HIGH) \
		$(RBEMAX2_PROTON_BLADDER_LOW) $(RBEMAX2_PROTON_BLADDER_MODE) $(RBEMAX2_PROTON_BLADDER_HIGH) \
		> log_$@.txt 2>&1


proton_bladder_sample_multiA_RBE_data.mat: merge_samples_multiA_RBE $(PROTON_BLADDER_FILELIST)
	./merge_samples_multiA_RBE $@ $(PROTON_BLADDER_FILELIST) && \
		rm -f $@.backup $(PROTON_BLADDER_FILELIST) $(foreach N,$(PROTON_BLADDER_FILELIST),log_$(N).txt)


$(PROTON_RECTUM_FILELIST): make_samples_multiA_RBE
	./make_samples_multiA_RBE $@ $(N_SAMPLES) 'data/VMATdvh/vmat%d.mat' 25 1 \
		'data/IMPTdvh/impt%d.mat' 25 1/1.1 Rectum \
		$(ALPHA_PROTON_RECTUM_MEAN) $(ALPHA_PROTON_RECTUM_STD) \
		$(BETA_PROTON_RECTUM_MEAN) $(BETA_PROTON_RECTUM_STD) \
		$(RBEMIN1_PROTON_RECTUM_LOW) $(RBEMIN1_PROTON_RECTUM_MODE) $(RBEMIN1_PROTON_RECTUM_HIGH) \
		$(RBEMAX1_PROTON_RECTUM_LOW) $(RBEMAX1_PROTON_RECTUM_MODE) $(RBEMAX1_PROTON_RECTUM_HIGH) \
		$(RBEMIN2_PROTON_RECTUM_LOW) $(RBEMIN2_PROTON_RECTUM_MODE) $(RBEMIN2_PROTON_RECTUM_HIGH) \
		$(RBEMAX2_PROTON_RECTUM_LOW) $(RBEMAX2_PROTON_RECTUM_MODE) $(RBEMAX2_PROTON_RECTUM_HIGH) \
		> log_$@.txt 2>&1


proton_rectum_sample_multiA_RBE_data.mat: merge_samples_multiA_RBE $(PROTON_RECTUM_FILELIST)
	./merge_samples_multiA_RBE $@ $(PROTON_RECTUM_FILELIST) && \
		rm -f $@.backup $(PROTON_RECTUM_FILELIST) $(foreach N,$(PROTON_RECTUM_FILELIST),log_$(N).txt)

