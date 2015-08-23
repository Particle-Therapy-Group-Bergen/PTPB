function samples = bootStrap(data, max_samples, sample_mode)
%samples = bootStrap(data [, max_samples, sample_mode])
%
% Calculates samples of the data vector by resampled it using the boot strap method.
%
% Parameters:
%   data - The data points to sample. Must be a vector.
%   max_samples - The maximum number of samples to generate. This is optional with
%                 a default value of 6435. If the data vector is small the function
%                 may produce less than max_samples number of samples, since all
%                 possible combinations will already have been resampled.
%   sample_mode - This is optional, but must be one of the following strings:
%                 'exhaustive' - Uses an exhaustive sampling method. max_samples
%                                will be ignored in this case. Note: Very slow for
%                                large data vectors.
%                 'random' - Uses a random sampling method to produce the samples.
%                            Will always produce exactly max_samples number of
%                            samples.
%                 'adaptive' - Will use the exhaustive method if the data vector is
%                              small, otherwise the random sampling method is used.
%
% The results are returned as a matrix of resampled vectors of the data vector.
% If data is a vector with length N then the resultant samples matrix will be MxN,
% with each sample placed in a matrix row.

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
    help bootStrap;
    return;
end

if ~ exist('max_samples')
  max_samples = 6435;
end
if ~ exist('sample_mode')
  sample_mode = 'adaptive';
end

switch sample_mode
  case 'exhaustive'
    samples = exhaustiveResample(data);
  case 'random'
    samples = randomResample(data, max_samples);
  case 'adaptive'
    % Check if we will produce no more than max_samples number of samples with
    % the exhaustive method. If so, then use that. Otherwise the random method.
    N = length(data);
    if factorial(2*N-1)/(factorial(N)*factorial(N-1)) <= max_samples
      samples = exhaustiveResample(data);
    else
      samples = randomResample(data, max_samples);
    end
  otherwise
    error('Unsupported sample_mode "%s".', sample_mode);
end
return;


function samples = randomResample(data, M)
% Perform random resampling of the data vector to produce M new resampled vectors.
N = length(data);
samples = zeros(M, N);
for m = 1:M
  samples(m,:) = data( round(rand(1, N) * N + 0.5) );
end
return;


function samples = exhaustiveResample(data)
% Performs exhaustive resampling of the data vector, i.e. finds every possible
% combination with replacement of data.
N = length(data);
resampleMatrix = combsrep(1:N, N);
samples = data(resampleMatrix);
return;
