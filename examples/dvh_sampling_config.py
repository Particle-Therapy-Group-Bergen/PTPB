# This is an example of a configuration file used with the sampleDVH.py tool.

# An optional list of organs. If no list of organs is given then all organs
# found in the input DVH files will be used.
#organs = ['Rectum']

# The following is an optional organ name mapping variable. It must be a
# dictionary mapping the name in DVH files (key) to the required name (value).
#organ_name_map = {'Bladder_P': 'Bladder', 'Rectum_P_MT': 'Rectum'}

# Optional volume bins at which to resample the DVHs.
#volume_bins = linspace(0,1,100)

# The interpolation method used with the interp1() function when resampling
# the DVHs to the volume-ratio bins.
interpolation_method = 'pchip'

# The dose bin uncertainty model for the DVH data points.
dose_binning_uncertainty = Box(0.05)

# The volume ratio bin uncertainty model for the DVH data points.
volume_ratio_uncertainty = Box(1e-6)

# Indicates the maximum number of samples to produce when random sampling in
# the boolStrap function. Refer to that function for more information.
bootstrap_max_samples = 6435

# The boot-strapping mode to use in the boolStrap function:
#   'adaptive', 'exhaustive' or 'random'
# Refer to that function for more information.
bootstrap_sample_mode = 'adaptive'
