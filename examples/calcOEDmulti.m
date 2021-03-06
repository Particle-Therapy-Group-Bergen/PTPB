function result = calcOEDmulti(filename, organ, integration_method, integration_tolerance)
%function result = calcOEDmulti(filename, organ, integration_method, integration_tolerance)
%
% This script is an example for calculating OED values for different organs found
% in a DVH file converted into a .mat workspace file. The converted file should contain
% the DVH_data variable.
% One can convert DVH files to this format using the convertDVHtoMatlabFile.sh script
% found in the DVHtools directory.
%
% This function takes the following parameters:
%  filename - The name of the .mat file to load and calculate from.
%  organ    - If a organ name is given then the calculation is only done for that organ.
%  integration_method - The integration method to use, if none is given then 'quad' is used.
%                       Refer to the ODE function for details about possible options.
%  integration_tolerance - The tollerance threshold to pass to the integration routine.
%
% Example usage to calculate for all organs using 'trapz' method and 1e-4 tollerance:
%
%  x = calcOEDmulti('data.mat', {}, 'trapz', 1e-4);
%
% or just for lungs:
%
%  x = calcOEDmulti('data.mat', 'Lungs', 'trapz', 1e-4);
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    Particle Therapy Project Bergen (PTPB) - tools and models for research in
%    cancer therapy using particle beams.
%
%    Copyright (C) 2014-2015 Particle Therapy Group Bergen
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

% Authors: Artur Szostak <artursz@iafrica.com>, Camilla H Stokkevaag <camilla.stokkevag@ift.uib.no>


% The following table remaps the names in the DVH files to standard ones:
%
%               Name in file          Name to map to
organnames = {
                'BODY',               'Body';

  };

% The following table maps different parameters to use for different model calculations.
%
%            The name of the    Threshold        Alpha        Competition model parameters
%            organ used in      parameter for    parameter                                         

%            the DVH files.     PlateauHall.     for LinExp.  alpha1  beta1       alpha2  beta2   (n fractions)	  delta
organtable = {
              {'Bladder_P',           4.5,          1.592,   0.00328,  0.00328/7.5,  0.25,  0.25/7.5,  28,	5.10 },
              {'Rectum_P_MT',         4.5,          0.240,   0.00984,  0.00984/5.4,  0.25,  0.25/5.4,  28,	0.26},

  };


% Build structures to map parameters more easily in the subsequent code.
organmap = {};
for n = 1:length(organtable)
  organmap{n*2-1} = organtable{n}{1};
  organmap{n*2  } = struct('threshold', organtable{n}{2}, 'alpha', organtable{n}{3},
                           'alpha1', organtable{n}{4}, 'beta1', organtable{n}{5},
                           'alpha2', organtable{n}{6}, 'beta2', organtable{n}{7},
                           'integrations', organtable{n}{8},'delta', organtable{n}{9});
end
organmap = struct(organmap{:});


% Load the data from file.
if ~ exist('filename')
  error('No data file name given.');
end
data = load(filename, 'DVH_data');


% Setup default for parameters that were not given.
if ~ exist('integration_method')
    integration_method = 'trapz';
end
if ~ exist('integration_tolerance')
    integration_tolerance = 1e-5;
end


% Process the organ structures:
doseResults = {};
resultCount = 1;
organs = data.DVH_data.structures;
one_organ_found = 0;
for n = 1:length(organs)
  s = organs{n};
  % Try map the name to a standard one.
  for m = 1:size(organnames, 1);
    if strcmp(s.structName, organnames{m,1})
      s.structName = organnames{m,2};
      break;
    end
  end
  if exist('organ')
    if ~ strcmp(organ, s.structName)
      continue;
    end
  end
  one_organ_found = 1;

  % Check if we can handle this organ.
  if ~ isfield(organmap, s.structName)
    continue;
  end
  if ~ isfield(s, 'dose')
    warning('The dose field could not be found so %s will be skipped.', s.structName);
    continue;
  end
  if ~ isfield(s, 'ratioToTotalVolume')
    warning('The ratioToTotalVolume field could not be found so %s will be skipped.', s.structName);
    continue;
  end
  if length(s.dose) == 0
    warning('The dose field is empty so %s will be skipped.', s.structName);
    continue;
  end
  if length(s.ratioToTotalVolume) == 0
    %NOTE: hacking the data to recover ratioToTotalVolume. Not clear if this procedure is correct!
    if ~ isfield(s, 'structureVolume')
      warning('The ratioToTotalVolume field is empty so %s will be skipped.', s.structName);
      continue;
    end
    warning('The ratioToTotalVolume field is empty so will use structureVolume for the %s instead.', s.structName);
    extraError = abs(s.structureVolume(1) - s.volume) / s.volume;
    s.ratioToTotalVolume = s.structureVolume / s.structureVolume(1);
  end
  if length(size(s.dose)) ~= length(size(s.ratioToTotalVolume))
    warning('The dose and ratioToTotalVolume fields are not the same size so %s will be skipped.', s.structName);
    continue;
  end
  if size(s.dose) ~= size(s.ratioToTotalVolume)
    warning('The dose and ratioToTotalVolume fields are not the same size so %s will be skipped.', s.structName);
    continue;
  end

  printf('Processing %s\n', s.structName);

  params = getfield(organmap, s.structName);
  threshold = params.threshold;
  alpha = params.alpha;
  alpha1 = params.alpha1;
  beta1 = params.beta1;
  alpha2 = params.alpha2;
  beta2 = params.beta2;
  integrations = params.integrations;
  delta = params.delta;

  % Rescale the volume fraction to ratio if set as percent.
  if max(s.ratioToTotalVolume) > 1
    s.ratioToTotalVolume = s.ratioToTotalVolume / 100;
  end
  datapoints = [s.dose; s.ratioToTotalVolume];

  % Here we perform the calculations and fill in the results. We try both linear
  % and pchip interpolation methods and vary the tolerance to get and estimate of
  % the numerical uncertainty.
  opts{1} = struct('integration_method', integration_method, 'tolerance', integration_tolerance, 'interpolation_method', 'pchip');
  opts{2} = struct('integration_method', integration_method, 'tolerance', integration_tolerance, 'interpolation_method', 'linear');
  opts{3} = struct('integration_method', integration_method, 'tolerance', integration_tolerance*10, 'interpolation_method', 'pchip');
  opts{4} = struct('integration_method', integration_method, 'tolerance', integration_tolerance*10, 'interpolation_method', 'linear');
  responseModels = {'Competition','CompNoF', 'LinPlat', 'PlateauHall', 'LNT'};
  modelResults = {};
  for k = 1:length(responseModels)
    responseModel = responseModels{k};
    switch responseModel
      case 'LNT'
        oed = @(opts) OED(responseModel, datapoints, opts);
      case 'PlateauHall'
        oed = @(opts) OED(responseModel, datapoints, opts, threshold);
      case 'LinExp'
        oed = @(opts) OED(responseModel, datapoints, opts, alpha);
      case 'Competition'
        oed = @(opts) OED(responseModel, datapoints, opts, alpha1, beta1, alpha2, beta2, n);
      case 'CompNoF'
        oed = @(opts) OED('LinExp', datapoints, opts, alpha2);
      case 'LinPlat'
        oed = @(opts) OED(responseModel, datapoints, opts, delta);
      otherwise
        error('Unsupported response model "%s".', responseModel);
    end
    doses = cellfun(oed, opts);  % apply oed to each set of options.
    dose = doses(1);   % select the best estimate for the dose.
    doseUncertainty = std(doses);  % estimate the uncertainty in the dose.
    if exist('extraError')
      doseUncertainty = sqrt(doseUncertainty^2 + (dose*extraError)^2);
    end
    printf('Calculated dose for %s = %g +/- %g\n', responseModel, dose, doseUncertainty);
    fflush(stdout);

    modelResults{k} = struct('dose', dose, 'doseUncertainty', doseUncertainty);
  end
  resultsRow = {};
  for k = 1:length(responseModels)
    resultsRow{2*k-1} = responseModels{k};
    resultsRow{2*k} = modelResults{k};
  end
  doseResults{resultCount} = s.structName;
  resultCount = resultCount + 1;
  doseResults{resultCount} = struct(resultsRow{:});
  resultCount = resultCount + 1;
end

if one_organ_found == 0
  printf('No organs found.\n');
end

result = struct(doseResults{:});
