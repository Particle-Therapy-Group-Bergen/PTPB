# This script is used to produce the mean DVH density plots.
# Note: the IMPT DVH's have to be scaled by dividing by 1.1.

BLADDER_ORGAN_ID="Bladder_P"
RECTUM_ORGAN_ID="Rectum_P_MT"
N_EXTRA_CURVES=10
QUANTILES="0.025:0.05:0.975"
NSAMPLES=100
BINNING="0:0.01:1"

if ! test -d data ; then
    echo "Please symlink the data directory to the local working path."
    echo "e.g.   ln -s <data-path> data"
    exit 1
fi
rm -rf ScaledIMPT
mkdir ScaledIMPT
applyResponseDVH.py -m '@(x) x./1.1' -V '0:0.001:1' -D '0:0.1:70' data/IMPTdvh/impt*.mat ScaledIMPT/
sampleDVH.py --nsamples "$NSAMPLES" --bins "$BINNING" --organ "$BLADDER_ORGAN_ID" data/CionDataPhysicalDose/HUH*.mat -O samples_cion_bladder.mat
sampleDVH.py --nsamples "$NSAMPLES" --bins "$BINNING" --organ "$BLADDER_ORGAN_ID" ScaledIMPT/impt*.mat -O samples_impt_bladder.mat
sampleDVH.py --nsamples "$NSAMPLES" --bins "$BINNING" --organ "$BLADDER_ORGAN_ID" data/VMATdvh/vmat*.mat -O samples_vmat_bladder.mat
plotMeanDVH.py --density --gnuplot --quantiles "$QUANTILES" --extra-curves "$N_EXTRA_CURVES" --outfile bladder.png samples_cion_bladder.mat samples_vmat_bladder.mat samples_impt_bladder.mat
sampleDVH.py --nsamples "$NSAMPLES" --bins "$BINNING" --organ "$RECTUM_ORGAN_ID" data/CionDataPhysicalDose/HUH*.mat -O samples_cion_rectum.mat
sampleDVH.py --nsamples "$NSAMPLES" --bins "$BINNING" --organ "$RECTUM_ORGAN_ID" ScaledIMPT/impt*.mat -O samples_impt_rectum.mat
sampleDVH.py --nsamples "$NSAMPLES" --bins "$BINNING" --organ "$RECTUM_ORGAN_ID" data/VMATdvh/vmat*.mat -O samples_vmat_rectum.mat
plotMeanDVH.py --density --gnuplot --quantiles "$QUANTILES" --extra-curves "$N_EXTRA_CURVES" --outfile rectum.png samples_cion_rectum.mat samples_vmat_rectum.mat samples_impt_rectum.mat
