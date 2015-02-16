function y = doseInterpolate(x, datapoints, method)
%y = doseInterpolate(x, datapoints)
%y = doseInterpolate(x, datapoints, method)
%
% This function interpolates discrete data points from a cumulative dose distribution into a continuous distribution.
%
%Where,
% x is volume fraction in the range [0..1]. x can be any of scalar, vector or matrix.
%
% datapoints are the discrete data points from the cumulative dose distribution.
% i.e. the volume fraction (in the range [0..1] on y-axis) as a function of dose (in Gray on x-axis).
% datapoints should be a 2 by N or N by 2 matrix.
%
% method is an optional parameter indicating the interpolation method to use in interp1().
% If not given 'pchip' is used by default.
%
% The output y is the dose for a given input x and has the same dimensions as x (if x is vector, y becomes a vector...).
%
%Example:
% xp = 0:0.1:50;
% yp = (1+erf(10-xp))/2;
% data = [xp; yp];
% x = 0:0.01:1;
% y = doseInterpolate(x, data, 'pchip');
% plot(x, y);
% title('Dose distribution');
% xlabel('volume fraction [volume element / total volume]');
% ylabel('dose [Gy]');
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    Particle Therapy Project Bergen (PTPB) - tools and models for research in
%    cancer therapy using particle beams.
%
%    Copyright (C) 2013 Particle Therapy Group Bergen
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
    help doseInterpolate;
    return;
end

% checking dimension of datapoints and returns error if not 2xN or Nx2
dim = size(datapoints);
if length(dim) ~= 2 || ~(dim(1) == 2 || dim(2) == 2)
    error('Parameter datapoints must be a 2 by N or N by 2 matrix.');
end

% Note: we actually need to deal with the function: y = f(x)
% where y is the dose and x is the volume fraction. However, the input
% datapoints are for the inverse: x = g(y). Thus, need to swap around the axes.
if dim(1) == 2
    % handling 2xN case
    yp = datapoints(1,:);
    xp = datapoints(2,:);
else
    % handling Nx2 case
    yp = datapoints(:,1);
    xp = datapoints(:,2);
end

% We need to remove any data points from the ends of the distribution that lie
% on the same x-axis coordinate, otherwise the interp1 function will fail.
starti = 1;
while starti < length(xp)-1 && xp(starti+1) == xp(starti)
    starti++;
end
endi = length(xp);
while endi > 2 && xp(endi-1) == xp(endi)
    endi--;
end

xp = xp(starti:endi);
yp = yp(starti:endi);

% Handle special case where the graph is f(y) == 0. All points are eliminated
% in the previous step.
if length(xp) < 2
    y = zeros(size(x));
    return;
end

% Handle another special case where dose is likely zero.
if length(xp) == 2 && xp(1) == 1 && xp(2) == 0 && yp(1) == 0 && yp(2) >= 0
    y = zeros(size(x));
    return;
end

% Now handle the points in the middle of the range that have two different y
% values but the same x value. Will need to separate the points along the
% x-axis by a small amount.
% Note: we only check the interior points in the range 2..n-1 since the previous
% steps guarantee that the first and second data points are different and the
% last two are also different.
% The first step is to get rid of the middle points for 3 or more consecutive
% duplicate points.
index = [1];
for n = 2:length(xp)-1
    if xp(n-1) ~= xp(n) || xp(n) ~= xp(n+1)
        index = [index, n];
    end
end
index = [index, length(xp)];
xp = xp(index);
yp = yp(index);
% The second step is to shift the remaining consecutive points by a small x.
K = 1e-7;
for n = 2:length(xp)-2
    if xp(n) == xp(n+1)
        xp(n)   = xp(n-1)*K + (1-K)*xp(n);
        xp(n+1) = xp(n+2)*K + (1-K)*xp(n+1);
    end
end

if ~ exist('method')
    % Select the interpolation method if not already given.
    % Note: spline and cubic seem to give numerical instabilities on the edges
    % of the x range so we stick to 'pchip' by default which appears to behave
    % better.
    method = 'pchip';
end
extrap = 'extrap';
y = interp1(xp, yp, x, method, extrap);
