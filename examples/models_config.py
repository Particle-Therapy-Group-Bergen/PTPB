# This is and example of an organ/model configuration file to use with the
# processPatients.py tool.

# Optional list of organs. If no list of organs is given then all organs found
# in the input DVH files will be used.
#organs = ['Bladder', 'Lungs']

# Optional list of models. If no list of models is given then all response
# models for which parameters are provided will be calculated.
#models = ['LNT', 'LinExp']

# The following is an optional organ name mapping variable. It must be a
# dictionary mapping the name in DVH files (key) to the required name (value).
organ_name_map = {'Bladder_P': 'Bladder', 'Rectum_P_MT': 'Rectum'}

# List of integration methods used for the OED calculation. This is a dictionary
# of method names as keys and the values indicate integration tolerances.
# See the method names in Octave for more details about the tolerance handling.
integration_methods = {'trapz': 1e-5,
                       'quad': 1e-5,
                       'quadv': 1e-5,
                       'quadl': 1e-5,
                       'quadgk': 1e-5}

# The interpolation methods used with the interp1() function when calculating
# the OED values.
interpolation_methods = ['linear', 'pchip', 'cubic', 'spline']

# The uncertainty model for the PlateauHall response model calculation used in
# the OED.m function.
# Available uncertainty models are:
#   Delta(value)
#   DoubleDelta(a, b)
#   Box(min, max)
#   Box95(low, high)
#   Triangle(min, mode, max)
#   Triangle95(low, mean, high)
#   Triangle95mode(low, mode, high)
#   Gaus(mean, sigma)
#   Gaus95(low, high)
#   LogNorm(mu, sigma)
#   LogNorm95(low, high)
# Refer to the sampleDistribution.m function for details about these
# distributions.
plateau_threshold_uncertainty = Triangle(4.0, 4.5, 5.0)

# Linear-exponential model alpha parameters per organ:
# Taken from table 1, Schneider et. al. 2005 paper:
#   "Estimation of radiation-induced cancer from three-dimensional dose
#    distributions: Concept of organ equivalent dose"
#   International journal of radiation oncology, biology, physics,
#   1 April 2005 (volume 61 issue 5 Pages 1510-1515
#   DOI: 10.1016/j.ijrobp.2004.12.040)
linexp_alphas = [
#                              Uncertainty
#    Organ name    Value     Low       High
    ('Stomach',    0.149,   -0.013,   0.014),
    ('Colon',      0.240,   -0.061,   0.075),
    ('Rectum',     0.240,   -0.061,   0.075),
    ('Bone',       0.033,   -0.027,   0.043),
    ('Liver',      0.487,   -0.252,   0.385),
    ('Lungs',      0.129,   -0.012,   0.016),
    ('Bladder',    1.592,   -0.356,   0.434),
    ('Thyroid',    0.033,   -0.013,   0.014),
    ('Prostate',   0.804,   -0.496,   0.622),
# Breast parameters taken from Schneider 2011, no CI available.
    ('Breast',     0.041,   -0.041,   0.041)
]

# The following is necessary code to setup the OrganParams objects from the
# linexp_alphas table and plateau_threshold_uncertainty parameter.
for organ, alpha, dlow, dhigh in linexp_alphas:
    OrganParams(
            name = organ,
            LNT = [],
            PlateauHall = [plateau_threshold_uncertainty],
            LinExp = [Triangle95mode(alpha + dlow, alpha, alpha + dhigh)]
        )

# The dose bin uncertainty model for the DVH data points.
dose_binning_uncertainty = Box(0.05)

# The volume ratio bin uncertainty model for the DVH data points.
volume_ratio_uncertainty = Box(1e-6)
