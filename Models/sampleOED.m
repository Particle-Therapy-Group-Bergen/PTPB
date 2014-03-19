function result = sampleOED(nsamples, filename, organs, models, options)
%result = sampleOED(nsamples, filename, organs, models, options)
%
%Samples the uncertainty distributions of the dose volume histograms and model parameters
%and calculates organ equivalent doses (OEDs) for these sampled values.
%i.e. the data points are jittered around their nominal values before calculating the
%OED to estimate the final uncertainty on the full OED calculation.
%
%Parameters:
% nsamples - The number of OED samples to produce.
% filename - The name of the converted DVH (.mat) data file containing the dose volume histograms.
% organs - Optional list of organs for which to calculate OED values. If none is explicitly given
%       then OEDs are calculated for all organs for which model parameters are available.
% models - Optional response models to use for the OED calculation. If none are explicitly given
%       then all possible models are computed. The valid values for models are any of:
%          'LNT' - Linear response model with no threshold.
%          'PlateauHall' - Linear model with a flat plateau threshold.
%          'LinExp' - Linear-exponential model.
%          'Competition' - Competition model (extension of LinExp) with only alpha and alpha/beta
%                         ratio parameters specified.
% options - A structure containing various parameters to control the calculation.
%   The valid field names of this structure are:
%      'organ_name_map' - A Nx2 cell matrix of key-value pairs to map organ names found in
%                         the data files to standard names. The first column contains the keys
%                         and the second the values. (default none)
%      'print_progress' - Integer indicating if progress information should be printed.
%                         Higher values give more verbosity. (default 0)
%      'debug_function' - Boolean indicating if extra debug checking should be performed. (default 0)
%      'struct_output' - Boolean value to force output to be in a full structure even if only dealing
%                         with one organ or model. (default 0)
%      'dosevolume_uncertainty_model' - Uncertainty distribution model to use for the dose volume
%                         histogram data points. Possible values: 'box', 'triangle' (default 'box')
%      'parameter_uncertainty_model' - Uncertainty distribution model to use for the model parameters.
%                         Possible values: 'box', 'triangle', 'gaussian' (default 'box')
%      'integration_method' - The integration method to use. Refer to "help OED" for a list of valid
%                         options. (default 'trapz')
%      'integration_tolerance' - The error tolerance level to use for the integration. (default 1e-4)
%      'interpolation_method' - The interpolation method to use. Refer to "help interp1" for a list of
%                         valid values. (default 'linear')
%      'dose_fractions' - The dose fraction value to use for the Competition models.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    Particle Therapy Project Bergen (PTPB) - tools and models for research in
%    cancer therapy using particle beams.
%
%    Copyright (C) 2014 Particle Therapy Group Bergen
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Authors: Artur Szostak <artursz@iafrica.com>

if ~ exist('nsamples') || nsamples < 1
  error('Must provide the number of samples to generate with a value >= 1.');
end
if ~ exist('filename')
  error('A file name to a converted DVH file (.mat file) must be provided.');
end

% Unpack the options structure and set defaults if needed.
print_progress = 0;
debug_function = 0;
struct_output = 0;
dosevolume_uncertainty_model = 'box';
parameter_uncertainty_model = 'box';
integration_method = 'trapz';
integration_tolerance = 1e-4;
interpolation_method = 'linear';
dose_fractions = 1;
if exist('options')
  if isfield(options, 'organ_name_map')
    organ_name_map = options.organ_name_map;
  end
  if isfield(options, 'print_progress')
    print_progress = options.print_progress;
  end
  if isfield(options, 'debug_function')
    debug_function = options.debug_function;
  end
  if isfield(options, 'struct_output')
    struct_output = options.struct_output;
  end
  if isfield(options, 'dosevolume_uncertainty_model')
    dosevolume_uncertainty_model = options.dosevolume_uncertainty_model;
  end
  if isfield(options, 'parameter_uncertainty_model')
    parameter_uncertainty_model = options.parameter_uncertainty_model;
  end
  if isfield(options, 'integration_method')
    integration_method = options.integration_method;
  end
  if isfield(options, 'integration_tolerance')
    integration_tolerance = options.integration_tolerance;
  end
  if isfield(options, 'integration_method')
    integration_method = options.integration_method;
  end
  if isfield(options, 'dose_fractions')
    dose_fractions = options.dose_fractions;
  end
end

if exist('organs') && ischar(organs)
  organs = {organs};
end
if exist('models') && ischar(models)
    models = {models};
end
if ~ exist('models') || length(models) == 0
  models = {'LNT', 'PlateauHall', 'LinExp', 'Competition'};
end

% Setup options for OED function.
oed_opts = struct(
    'integration_method', integration_method,
    'tolerance', integration_tolerance,
    'interpolation_method', interpolation_method
  );

% Fetch the dose volume histograms, taking care of passing the right set of parameters.
func_params = {filename};
if exist('organ_name_map')
  func_params = {func_params{:}, organ_name_map};
end
if exist('organs')
  func_params = {func_params{:}, organs{:}};
end
doseVolumeHistos = getDoseVolumeHistogram(func_params{:});

if debug_function
  % Check that the ranges of the dose volume histogram data points are correct if we are debugging.
  names = fieldnames(doseVolumeHistos);
  for k = 1:length(names)
    dose = doseVolumeHistos.(names{k}).datapoints(:,1);
    volumeRatio = doseVolumeHistos.(names{k}).datapoints(:,2);
    doseMin = doseVolumeHistos.(names{k}).range_low(:,1);
    volumeRatioMin = doseVolumeHistos.(names{k}).range_low(:,2);
    doseMax = doseVolumeHistos.(names{k}).range_high(:,1);
    volumeRatioMax = doseVolumeHistos.(names{k}).range_high(:,2);
    for n = 1:length(dose)
      if ~ (doseMin(n) <= dose(n) && dose(n) <= doseMax(n))
        error('Order assertion failure: doseMin(%d) = %.16e, dose(%d) = %.16e, doseMax(%d) = %.16e',
              n, doseMin(n), n, s.dose(n), n, doseMax(n));
      end
      if ~ (volumeRatioMin(n) <= volumeRatio(n) && volumeRatio(n) <= volumeRatioMax(n))
        error('Order assertion failure: volumeRatioMin(%d) = %.16e, volumeRatio(%d) = %.16e, volumeRatioMax(%d) = %.16e',
              n, volumeRatioMin(n), n, volumeRatio(n), n, volumeRatioMax(n));
      end
    end
  end
end

% Fetch the model parameters for the different organs.
params = getParameters();

names = fieldnames(doseVolumeHistos);
for k = 1:length(names)
  % Check if we have any model parameters for the given organ at all.
  organ_name = names{k};
  if ~ isfield(params, organ_name)
    continue
  end
  mp = params.(organ_name);
  % Go through the different models and if we have parameters available for
  % the given organ and model combination then perform the computation.
  for n = 1:length(models)

    % First prepare a vector of model specific parameters and their ranges for the OED function.
    switch models{n}
      case 'LNT'
        mp_values = [];
        mp_range_low = [];
        mp_range_high = [];
        mp_constants = [];
      case 'PlateauHall'
        if isfield(mp, 'plateau_threshold')
          mp_values = [mp.plateau_threshold.value];
          mp_range_low = [mp.plateau_threshold.range_low];
          mp_range_high = [mp.plateau_threshold.range_high];
          mp_constants = [];
        else
          continue;
        end
      case 'LinExp'
        if isfield(mp, 'linexp_alpha')
          mp_values = [mp.linexp_alpha.value];
          mp_range_low = [mp.linexp_alpha.range_low];
          mp_range_high = [mp.linexp_alpha.range_high];
          mp_constants = [];
        else
          continue;
        end
      case 'Competition'
        if isfield(mp, 'competition_alpha1') && isfield(mp, 'competition_alpha2') && isfield(mp, 'competition_alpha_beta_ratio')
          mp_values = [mp.competition_alpha1.value, mp.competition_alpha2.value, mp.competition_alpha_beta_ratio.value];
          mp_range_low = [mp.competition_alpha1.range_low, mp.competition_alpha2.range_low, mp.competition_alpha_beta_ratio.range_low];
          mp_range_high = [mp.competition_alpha1.range_high, mp.competition_alpha2.range_high, mp.competition_alpha_beta_ratio.range_high];
          mp_constants = [dose_fractions];
        else
          continue;
        end
      otherwise
        error('Model "%s" is not supported.', models{n});
    end

    % Jitter the model parameters to generate nsamples worth of calculated doses.
    mp_values = repmat(mp_values, nsamples, 1);
    mp_range_low = repmat(mp_range_low, nsamples, 1);
    mp_range_high = repmat(mp_range_high, nsamples, 1);
    mp_jittered_values = sampleParameters(mp_values, mp_range_low, mp_range_high, parameter_uncertainty_model);

    % Now produce the required number of sampled OED points.
    dose = [];
    for m = 1:nsamples

      % Jitter the histogram data points.
      datapoints = sampleDoseVolumeHisto(doseVolumeHistos.(names{k}).datapoints,
                                         doseVolumeHistos.(names{k}).range_low,
                                         doseVolumeHistos.(names{k}).range_high,
                                         dosevolume_uncertainty_model);
      if debug_function
        % Check the data points are in the correct order if we are debugging:
        for j = 2:size(datapoints, 1)
          p1 = datapoints(j-1,:);
          p2 = datapoints(j  ,:);
          x1 = p1(2); y1 = p1(1);
          x2 = p2(2); y2 = p2(1);
          if x2 > x1 || y1 > y2
            error('Data point order check failed: n = %d, x1 = %.20e, x2 = %.20e, y1 = %.20e, y2 = %.20e\n', n, x1, x2, y1, y2);
          end
        end
        % Plot the points for debugging:
        oldpoints = doseVolumeHistos.(names{k}).datapoints;
        plot(oldpoints(:,2), oldpoints(:,1), '+', datapoints(:,2), datapoints(:,1), 'o');
      end

      if isempty(mp_jittered_values)
        model_params = {};
      else
        model_params = [mp_jittered_values(m,:), mp_constants];
        % Convert model parameter matrix to a cell array to pass to the calcOED function.
        model_params = mat2cell(model_params, 1, ones(length(model_params), 1));
      end

      dose(m) = calcOED(models{n}, datapoints, oed_opts, model_params{:});

      if print_progress > 2
        printf('Dose = %f.\n', dose(m));
        fflush(stdout);
      end
    end

    % Add dose vector to results.
    result.(organ_name).(models{n}) = dose';

    if print_progress > 1
      printf('Completed model "%s".\n', models{n});
      fflush(stdout);
    end
  end

  if print_progress > 0
    printf('Completed organ "%s".\n', organ_name);
    fflush(stdout);
  end
end

if struct_output == 0
  % Return simplified format if only one field exists in the structure.
  if length(models) == 1
    names = fieldnames(result);
    for n = 1:length(names)
      subresults = result.(names{n});
      result.(names{n}) = subresults.(fieldnames(subresults){1});
    end
  end
  if length(fieldnames(result)) == 1
    result = result.(fieldnames(result){1});
  end
end

return;


function dose = calcOED(responseModel, datapoints, options, varargin)
% Calculates the organ equivalent dose via the OED function, making all necessary pre-calculations.
switch responseModel
  case 'Competition'
    alpha1 = varargin{1};
    alpha2 = varargin{2};
    beta1 = alpha1 / varargin{3};
    beta2 = alpha2 / varargin{3};
    n = varargin{4};
    dose = OED(responseModel, datapoints, options, alpha1, beta1, alpha2, beta2, n);
  otherwise
    dose = OED(responseModel, datapoints, options, varargin{:});
end
return;


function result = sampleDoseVolumeHisto(datapoints, range_low, range_high, uncertainty_model)
% Samples the uncertainty distributions for the dose volume histogram data points.
% i.e. takes the data points and jitters them by a random error taken from a probability
% distribution given by the specified model.

switch uncertainty_model
  case 'triangle'
    result = randTriangle(range_low, range_high, datapoints);
  case 'box'
    result = randBox(range_low, range_high);
  otherwise
    error('The uncertainty distribution model "%s" is not supported.', uncertainty_model);
end
return;


function result = sampleParameters(value, range_low, range_high, uncertainty_model)
% Samples the uncertainty distributions for a parameter set, similar to sampleDoseVolumeHisto.

switch uncertainty_model
  case 'triangle'
    result = randTriangle(range_low, range_high, value);
  case 'box'
    result = randBox(range_low, range_high);
  case 'gaussian'
    result = randn(size(value)) .* 0.5 .* abs(range_high - range_low) + value;
  otherwise
    error('The uncertainty distribution model "%s" is not supported.', uncertainty_model);
end
return;


function X = randTriangle(a, b, c)
% Generate random number matrix with the same dimentions as a, distributed according to
% a triangular distribution with min a, max b and mode (peak) c.

U = rand(size(a));
F = (b-a).*U < (c-a);
X = (a + sqrt(U.*(b-a).*(c-a))) .* F + (b - sqrt((1-U).*(b-a).*(b-c))).*(~F);
return;


function X = randBox(a, b)
% Generate random number matrix with the same dimentions as a, distributed according to
% a box distribution with min a, max b.

U = rand(size(a));
X = (b-a).*U + a;
return;
