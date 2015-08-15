function dVdD = doseDifferencial(dose, doseCumulative, interpMethod, tolerance)
%dVdD = doseDifferencial(dose, doseCumulative [, interpMethod, tolerance])
%
% Calculates the differencial of the cumulative dose volume histogram, which must
% be loaded from file with the getDoseVolumeHistogram() function.
%
%Where,
% dose is the vector of dose values at which to calculate the differencial
% function. i.e. the X values of the DVH function (with dose on the X-axis).
%
% doseCumulative are the cumulative dose distribution data points as a function
% of dose (dose on X-axis). This should be the two column matrix as loaded by
% getDoseVolumeHistogram() function.
%
% interpMethod indicates the interpolation method to use for the interp1 function.
% by default this is 'pchip'.
%
% tolerance is the numerical tolerance value used in the numerical derivative
% calculation. i.e. h in the following equation:
%    df   f(x+h) - f(x-h)
%    -- = ---------------
%    dx         2h
% By default this value is set to 1e-9.
%
% The output is a vector of volume fraction differencial values.
%
%Example:
% organ = 'Bladder';
% dvh = getDoseVolumeHistogram('mydvhfile.mat', organ);
% data = dvh.(organ).datapoints;
% x = 0:0.1:50;
% y = doseDifferencial(x, data, 'pchip');
% plot(x, y);
% title('Dose differencial');
% xlabel('dV/dD');
% ylabel('dose [Gy]');
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    Particle Therapy Project Bergen (PTPB) - tools and models for research in
%    cancer therapy using particle beams.
%
%    Copyright (C) 2013-2014 Particle Therapy Group Bergen
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

if nargin == 0
    % Print help message if no arguments are given.
    help doseDifferencial;
    return;
end

if ~ exist('interpMethod')
  interpMethod = 'pchip';
end
if ~ exist('tolerance')
  tolerance = 1e-9;
end

y2 = interpolate(doseCumulative(:,1), doseCumulative(:,2), dose+tolerance, interpMethod);
y1 = interpolate(doseCumulative(:,1), doseCumulative(:,2), dose-tolerance, interpMethod);
dVdD = (y1-y2) / (2*tolerance);
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
