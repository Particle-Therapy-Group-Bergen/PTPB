
PATIENT = 12
TOLERANCE = 1e-4
ALPHA = 0.25
BETA = 0.033
RBE_MIN = 1.25
RBE_MAX = 6
OUTFILE = sample_data.mat
ITERATIONS = 32
N_SAMPLES = 64

OCTAVE = octave -q --path ~/bin/PTPB_mfiles


FILEIDS = $(shell N=1; while test $$N -le $(ITERATIONS); do echo $$N ; N=$$((N+1)) ; done)
FILELIST = $(foreach N,$(FILEIDS),samples$(N).mat)

.PHONY: all clean merge

all: merge

clean:
	rm -rf make_samples merge_samples $(FILELIST)


define merge_samples
#!/bin/sh
outfile="$$1"
shift
files=""
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
export merge_samples

merge_samples: Makefile
	echo "$$merge_samples" > $@
	chmod +x $@


merge: merge_samples $(FILELIST)
	if test -f $(OUTFILE) ; then \
		mv $(OUTFILE) $(OUTFILE).backup; \
		./merge_samples $(OUTFILE) $(OUTFILE).backup $(FILELIST) \
			&& rm -f $(OUTFILE).backup $(FILELIST); \
	else \
		./merge_samples $(OUTFILE) $(FILELIST) \
			&& rm -f $(FILELIST); \
	fi


define make_samples
#!/bin/sh
exec $(OCTAVE) <<EOF
Nsamples = $$2;
patients = [$$3];
filepat1 = 'data/VMATdvh/vmat%d.mat';
n1 = 25;
scale1 = 1;
filepat2 = 'data/CionDataPhysicalDose/HUH%dphysical_dvh.mat';
n2 = 12;
scale2 = 1;
organ = 'Bladder';
alpha_distrib = struct('type', 'delta', 'params', {{$$4}});
beta_distrib = struct('type', 'delta', 'params', {{$$5}});
RBEmin_distrib = struct('type', 'delta', 'params', {{$$6}});
RBEmax_distrib = struct('type', 'delta', 'params', {{$$7}});
opts = struct('integration_method', 'trapz',
              'integration_tolerance', $$8,
              'interpolation_method', 'pchip',
              'sample_dvh', 1,
              'bootstrap_samples', 1,
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
export make_samples

make_samples: Makefile
	echo "$$make_samples" > $@
	chmod +x $@

$(FILELIST): make_samples
	./make_samples $@ $(N_SAMPLES) $(PATIENT) $(ALPHA) $(BETA) $(RBE_MIN) $(RBE_MAX) $(TOLERANCE)

