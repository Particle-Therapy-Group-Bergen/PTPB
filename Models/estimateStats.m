function result = estimateStats(x, nsamples)
%result = estimateStats(x, nsamples)
%
% Calculates the minimum, maximum, mean, median and standard deviation of a
% distribution using the bootstrapping method. i.e resampling the distribution
% with replacement for 'nsamples' number of times to get a distribution of
% the above mentioned statistics. From these re-samples we can estimate the
% value of the statistics themselves and their uncertainties.
%
% Parameters:
%   x - The data points of the distribution.
%   nsamples - The number of times to resample the data points x. The more
%              times one resamples the distribution the more precision is
%              achieved in the final results.
%
% The results are returned as a structure with each field representing one
% of the statistics estimated for the distribution x.

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

xlen = length(x);

samples = [];
for n = 1:nsamples
  y = x( round(rand(1, xlen) * xlen + 0.5) );
  samples(n,1:5) = [min(y), max(y), mean(y), median(y), std(y)];
end

ymean = mean(samples);
ysigma = std(samples) ./ sqrt(nsamples);
field = {'min', 'max', 'mean', 'median', 'stdev'};

result = struct;
for n = 1:5
  result.(field{n}) = struct('value', ymean(n), 'uncertainty', ysigma(n));
end
