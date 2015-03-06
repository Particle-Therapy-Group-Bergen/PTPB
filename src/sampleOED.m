function result = sampleOED(dvh, N, responseModel, integrationMethods, interpolationMethods, varargin)
%result = sampleOED(dvh, N, responseModel, integrationMethods, interpolationMethods, ...)
%
% This function randomly samples input dose volume histogram and response model
% parameter uncertainty distributions to calculate OED values. This monte-carlo
% method produces an output distribution of OED samples from which statistical
% values can be calculated, such as the likely OED value and confidence
% intervals.
%
%Parameters:
%
% dvh - The dose volume histogram structure as returned by the
%       getDoseVolumeHistogram function. With two additional fields:
%       dose_binning_uncertainty_model and volume_ratio_uncertainty_model which
%       should be the uncertainty models passed to sampleDistribution().
%
% N - The number of samples to produce.
%
% responseModel - A string indicating the response model to use. Refer to the
%                 OED function for details.
%
% integrationMethods - A cell array of strings naming the integration methods to
%                      use for the OED function. Will randomly select between
%                      methods before calculating each new OED sample.
%
% interpolationMethods - A cell array of strings naming the interpolation
%                        methods to use for the OED function. Will randomly
%                        select between methods before calculating each new OED
%                        sample.
%
% Additional parameters should be structures, one for each response model
% parameter that must be passed to the OED function. Each structure contains
% the following two fields to define the uncertainty distribution of the
% parameter:
%   uncertainty_model - A string indicating a distribution model to use. Can be
%                       any value that is accepted by the function
%                       sampleDistribution.
%   params - A cell array of parameters passed to sampleDistribution that define
%            the distribution. The exact number of parameters and their meaning
%            depends on the value of uncertainty_model chosen. For details refer
%            to the sampleDistribution function.
%
% This function returns a vector of OED samples.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    Particle Therapy Project Bergen (PTPB) - tools and models for research in
%    cancer therapy using particle beams.
%
%    Copyright (C) 2015 Particle Therapy Group Bergen
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

% Author: Artur Szostak <artursz@iafrica.com>

if nargin == 0
    % Print help message if no arguments are given.
    help sampleOED;
    return;
end

% Sample the DVH N times.
dose = dvh.datapoints(:,1);
doseLow = dvh.range_low(:,1);
doseHigh = dvh.range_high(:,1);
volume = dvh.datapoints(:,2);
volumeLow = dvh.range_low(:,2);
volumeHigh = dvh.range_high(:,2);
switch dvh.dose_binning_uncertainty_model
    case 'triangle'
        dose_samples = sampleDistribution(dvh.dose_binning_uncertainty_model, N, doseLow', dose', doseHigh');
    case 'triangle95'
        dose_samples = sampleDistribution(dvh.dose_binning_uncertainty_model, N, doseLow', dose', doseHigh');
    otherwise
        dose_samples = sampleDistribution(dvh.dose_binning_uncertainty_model, N, doseLow', doseHigh');
end
switch dvh.volume_ratio_uncertainty_model
    case 'triangle'
        volume_samples = sampleDistribution(dvh.volume_ratio_uncertainty_model, N, volumeLow', volume', volumeHigh');
    case 'triangle95'
        volume_samples = sampleDistribution(dvh.volume_ratio_uncertainty_model, N, volumeLow', volume', volumeHigh');
    otherwise
        volume_samples = sampleDistribution(dvh.volume_ratio_uncertainty_model, N, volumeLow', volumeHigh');
end

% Sample the response model parameters N times.
if length(varargin) > 0
    modelParams = [];
    for n = 1:length(varargin)
        p = varargin{n};
        modelParams(:,n) = sampleDistribution(p.uncertainty_model, N, p.params{:});
    end
else
    modelParams = zeros(N, 0);
end

result = zeros(N, 1);
for n = 1:N
    dvh_samples = [dose_samples(n,:); volume_samples(n,:)]';
    [nr, nc] = size(modelParams);
    if nr*nc > 0
        model_param_samples = mat2cell(modelParams(n,:), 1, ones(length(varargin), 1));
    else
        model_param_samples = {};
    end

    % Randomly select the integration and interpolation method.
    index = randInt(length(integrationMethods));
    integrationMethod = integrationMethods{index};
    index = randInt(length(interpolationMethods));
    interpolationMethod = interpolationMethods{index};

    options = struct('integration_method', integrationMethod.name,
                     'tolerance', integrationMethod.tolerance,
                     'interpolation_method', interpolationMethod);

    result(n) = OED(responseModel, dvh_samples, options, model_param_samples{:});
end


function y = randInt(imax)
if exist('randi')
    y = randi(imax);
else
    y = floor(rand * imax) + 1;
end
return;
