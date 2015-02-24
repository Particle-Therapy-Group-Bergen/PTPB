function [dose, volratio] = calcDVHquantiles(samples, V, D, Q, interpMethod)
%[dose, volratio] = calcDVHquantiles(samples, V, D, Q [, interpMethod])
%
% This function calculates the quantiles for the DVH samples produced by the
% sampleDVH() function. The calculation is done along each volume-ratio bin.
%
%Parameters:
%
% samples - The matrix of sample values returned by the sampleDVH() function.
%
% V - The volume ratio bins (interpolation points) that were used in sampleDVH().
%
% D - The dose bins to interpolate the quantiles to.
%
% Q - A vector of quantile values to calculate.
%
% interpMethod - Optional interpolation method to use, default = 'pchip'.

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
    help calcDVHquantiles;
    return;
end

if ~ exist('interpMethod')
    % Select the interpolation method if not already given.
    interpMethod = 'pchip';
end

[nr, nc] = size(samples);
if length(V) ~= nc
    error('The number of columns in "samples" and the length of "V" must match.');
end

[nr, nc] = size(D);
if nr > nc
    D = D';
end
dose = repmat(D, length(Q), 1)';
q = quantile(samples, Q, 1, 8);
[nr, nc] = size(q);
qrows = mat2cell(q, ones(nr, 1), nc);
func = @(x) interpolate(x, V, D, interpMethod);
vols = cellfun(func, qrows, 'UniformOutput', false);
volratio = cell2mat(vols)';
return;


function y = interpolate(xp, yp, x, method)
% Interpolation function that sets all values below the data points X range to 1 and
% all above to zero.
y = interp1(xp, yp, x, method);
min_non_nan = min(find(~ isnan(y)));
indices = find(isnan(y));
if length(min_non_nan) == 1
    nans_below_range = indices(find(indices < min_non_nan));
    nans_above_range = indices(find(indices > min_non_nan));
    y(nans_below_range) = 1;
    y(nans_above_range) = 0;
else
    y = zeros(size(y));
end
return;
