function result = sampleLAR(oed_samples, age, interpolation_methods, varargin)
%result = sampleLAR(oed_samples, age, interpolation_methods, point1 [, point2, ...])
%
% This function samples coefficient LAR data points at a particular age and computes
% the LAR value given a set of OED samples.
%
%Parameters:
%
% oed_samples - The samples to multiply the sampled LAR coefficient with. These
%               should be created with sampleOED.m
%
% age - The age at which to calculate the LAR value.
%
% interpolation_methods - A cell array of strings naming the interpolation
%                         methods to use for the interp1 function. Will randomly
%                         select between methods before calculating each new LAR
%                         sample.
%
% point1, point2 ... - All additional function parameters are treated as LAR
%                      data points. These must each be 2 element cell arrays.
%                      The first cell array entry must give the age value of
%                      the data point and the second entry must be a structure
%                      indicating the uncertainty model for the LAR coefficient
%                      at that age point. The structure should have two fields
%                        'uncertainty_model' - The uncertainty model name.
%                                              Should be a value compatible
%                                              with sampleDistribution().
%                        'params' - A cell array of parameters for the given
%                                   uncertainty model, as passed to the
%                                   sampleDistribution() function.
%
% This function returns a vector of LAR samples.

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
    help sampleLAR;
    return;
end

% Construct a vector of ages from the LAR data points and a matrix of
% samples from the LAR uncertainty distributions. The samples will be aligned
% along the columns for each age.
ages = [];
lar_points = [];
N = length(oed_samples);
for n = 1:length(varargin)
    ages(n) = varargin{n}{1};
    p = varargin{n}{2};
    lar_points(:,n) = sampleDistribution(p.uncertainty_model, N, p.params{:});
end
[ages, I] = sort(ages);
lar_points = lar_points(:,I);

lar_samples = [];
for n = 1:N
    % Randomly select the interpolation method.
    index = randInt(length(interpolation_methods));
    interpolation_method = interpolation_methods{index};

    % Interpolate the LAR data points for the patient age.
    lar_samples(n) = interp1(ages, lar_points(n,:), age, interpolation_method);
end
if size(lar_samples) ~= size(oed_samples)
    result = lar_samples' .* oed_samples;
else
    result = lar_samples .* oed_samples;
end
return;


function y = randInt(imax)
if exist('randi')
    y = randi(imax);
else
    y = floor(rand * imax) + 1;
end
return;
