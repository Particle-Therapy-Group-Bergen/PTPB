function y = doseInterpolate(x,datapoints)
%y = doseInterpolate(x,datapoints)
%
% This function interpolates input discrete data points into a continuous distribution.
%
%Where,
% x is volume fraction in the range [0..1]. x can be any of scalar, vector or matrix.
%
% datapoints are the discrete data points from the cumulative dose distribution.
% i.e. the volume fraction (in the range [0..1] on y-axis) as a function of dose (in Gray on x-axis).
% datapoints should be a 2 by N or N by 2 matrix.
%
% The output y is the dose for a given input x and has the same dimensions as x (if x is vector, y becomes a vector...).

% checking dimension of datapoints and returns error if not 2xN or Nx2
dim=size(datapoints);
if length(dim) != 2 || !(dim(1)==2 || dim(2)==2)
	error('Parameter datapoints must be a 2 by N or N by 2 matrix.');
end

% Note: we actually need to deal with the function: y = f(x)
% where y is the dose and x is the volume fraction. However, the input datapoints are
% for the inverse: x = g(y). Thus, need to swap around the axes.
if dim(1)==2
	% handling 2xN case
	yp=datapoints(1,:);
	xp=datapoints(2,:);
else
	% handling Nx2 case
	yp=datapoints(:,1);
	xp=datapoints(:,2);
end

% We need to remove any data points from the ends of the distribution that lie on the
% same x-axis coordinate, otherwise the interp1 function will fail.
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

method='pchip';  % Note: spline and cubic seem to give numerical instabilities on the edges of the x range.
extrap='extrap';
y = interp1(xp, yp, x, method, extrap);
