function [result, samples] = estimateStats(x, nsamples, lower_quantile, upper_quantile)
%[result, samples] = estimateStats(x, nsamples, lower_quantile, upper_quantile)
%
% Calculates the minimum, maximum, lower quantile, upper quantile, mean, median
% and standard deviation of a distribution using the bootstrapping method. i.e
% resampling the distribution with replacement for 'nsamples' number of times
% and calculating the statistics on the new distribution. This gives us a
% number of samples for each statistic which itself forms a distribution. From
% these resampled distributions we can estimate the value of the statistics and
% their uncertainties.
% Note: by default the upper and lower quantile are chosen to form a 95%
% confidence interval.
%
% Parameters:
%   x - The data points of the distribution.
%   nsamples - The number of times to resample the data points x. The more
%              times one resamples the distribution the more precision is
%              achieved in the final results. (optional, default = 1000)
%   lower_quantile - The lower quantile value in the range [0..1].
%                    (optional, default = 0.025)
%   upper_quantile - The upper quantile value in the range [0..1].
%                    (optional, default = 0.975)
%
% The results are returned as a structure with each field representing one
% of the statistics estimated for the distribution x.
% In addition, the calculated samples are returned as a Nx7 matrix with the 7
% columns indicating, the minimum, maximum, lower quantile, upper quantile,
% mean, median and standard deviation in that order.

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

% Author: Artur Szostak <artursz@iafrica.com>

if nargin == 0
    % Print help message if no arguments are given.
    help estimateStats;
    return;
end

if ~ exist('nsamples')
  nsamples = 1000;
end
if ~ exist('lower_quantile')
  lower_quantile = 2.5e-2;
end
if ~ exist('upper_quantile')
  upper_quantile = 97.5e-2;
end

xlen = length(x);

samples = [];
for n = 1:nsamples
  y = x( round(rand(1, xlen) * xlen + 0.5) );
  q1 = quantile(y, lower_quantile, 1, 8);
  q2 = quantile(y, upper_quantile, 1, 8);
  samples(n,1:7) = [min(y), max(y), q1, q2, mean(y), median(y), std(y)];
end

ymean = mean(samples);
ysigma = std(samples) ./ sqrt(nsamples);
field = {'min', 'max', 'lower_quantile', 'upper_quantile', 'mean', 'median', 'stdev'};

result = struct;
for n = 1:length(field)
  result.(field{n}) = struct('value', ymean(n), 'uncertainty', ysigma(n));
end
