function effective_dvh = applyResponse(dvh, V, D, interpMethod, funcname, varargin)
%effective_dvh = applyResponse(dvh, V, D, interpMethod, funcname, varargin)
%
% Applies a response function calculation to a DVH and returns the effective
% DVH after the applied response.
%
%Parameters:
%
% dvh - The DVH structure as is returned by the getDoseVolumeHistogram() function.
%
% V - The points on which the volume ratio dimension should be interpolated.
%
% D - The points on which the dose dimension should be interpolated.
%
% interpMethod - The interpolation method to use for the interp1() function.
%
% funcname - The name of the dose response function to use or a function handle
%            to specify a response function directly.
%
% Any additional parameters will be passed onto the response function.

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
    % Print help message if no arguments were given.
    help applyResponse;
    return;
end

if ischar(funcname)
    % Select the response function based on the response model or assume it is
    % a name for an existing Matlab function.
    switch funcname
        case 'LNT'
            func = @(x) x;
        case 'PlateauHall'
            func = @(x) PlateauHall(x, varargin{:});
        case 'LinExp'
            func = @(x) LinExp(x, varargin{:});
        case 'Competition'
            func = @(x) Competition(x, varargin{:});
        case 'LinPlat'
            func = @(x) LinPlat(x, varargin{:});
        otherwise
            eval(sprintf('func = @(x) %s(x, varargin{:});', funcname));
    end
else
    func = @(x) funcname(x, varargin{:});
end

% Change vectors to column format.
[nr, nc] = size(D);
if nc > nr
    D = D';
end
[nr, nc] = size(V);
if nc > nr
    V = V';
end

% Interpolate the data points on bins along the volume ratio dimension. Next, apply
% the response function to the doses at these points. Finally, interpolate again
% on bins on the dose dimension to get back to a regular dose bin grid.
% Note we have to sort the calculated effective doses to get back a monotonic
% function since the response model may be non-monotonic.
doses = doseInterpolate(V, dvh.datapoints, interpMethod);
effectiveDoses = sort(func(doses), 'descend');
volumeRatios_datapoints = interpolate(effectiveDoses, V, D, interpMethod);
doses = doseInterpolate(V, dvh.range_low, interpMethod);
effectiveDoses = sort(func(doses), 'descend');
volumeRatios_range_low = interpolate(effectiveDoses, V, D, interpMethod);
doses = doseInterpolate(V, dvh.range_high, interpMethod);
effectiveDoses = sort(func(doses), 'descend');
volumeRatios_range_high = interpolate(effectiveDoses, V, D, interpMethod);

effective_dvh = struct('datapoints', [D volumeRatios_datapoints],
                       'range_low', [D volumeRatios_range_low],
                       'range_high', [D volumeRatios_range_high]);
return;


function y = PlateauHall(d, threshold)
y = d .* (d < threshold) + threshold .* (d >= threshold);
return;


function y = LinExp(d, alpha)
y = d.*exp(-alpha.*d);
return;


function y = Competition(d, alpha1, beta1, alpha2, beta2, n)
% n is the number of dose fractions.
y = (d + beta1./alpha1.*d.^2./n).*exp(-(alpha2.*d + beta2.*d.^2./n));
return;


function y = LinPlat(d, delta)
y = ((1-exp(-delta.*d))./delta);
return;


function y = interpolate(xp, yp, x, method)
% Interpolation function that sets all values below the data points X range to 1 and
% all above to zero.
% We first remove any duplicate X values since methods like pchip will fail otherwise.
xp2 = [xp(1)];
yp2 = [yp(1)];
k = 2;
for n = 2:length(xp)
    if xp(n) ~= xp(n-1)
        xp2(k) = xp(n);
        yp2(k) = yp(n);
        k += 1;
    end
end
xp = xp2;
yp = yp2;
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
