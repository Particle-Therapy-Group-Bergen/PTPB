function y = doseInterpolate(x,datapoints) 
% function interpolates input discrete datapoints into a continous distribution
% x is volume fraction, x can be any of scalar, vector or matrix
% datapoints are the discrete datapoints from the cummulative dose distribution
% datapoints should be a 2xn matrix
% y is the dose for a given x and it has the same dimesions as x (if x is vector, y becomes a vector...)

% checking dimension of datapoints and returns error if not 2xn or nx2  
dim=size(datapoints);
if length(dim) != 2 || !(dim(1)==2 || dim(2)==2)
	error('Datapoint must be a 2xn or nx2 matrix');
end

if dim(1)==2
	% handling 2xn case
	xp=datapoints(1,:);
	yp=datapoints(2,:);
else
	% handling nx2 case
	xp=datapoints(:,1);
	yp=datapoints(:,2);
end

% Calculate numerical derivative df/dx with: df/dx = (f(x+h) - f(x-h)) / (2*h)
% h can be adjusted as needed, controls accuracy and stability.
h = 1e-10;
method='spline';
fxph = interp1(xp, yp, x+h, method);  % fxph = f(x+h)
fxmh = interp1(xp, yp, x-h, method);  % fxmh = f(x-h)
y = (fxph - fxmh) ./ (2.*h);

