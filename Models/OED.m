function dose = OED(responseModel, doseCumulative, varargin)
%dose = OED(responseModel, doseCumulative, ...)
%
% Calculates the Organ Equivalent Dose (OED) from the cumulative dose distribution data points.
%
%Where,
% responseModel can be either 'LNT', 'PlateauHall' or 'LinExp'.
%
% doseCumulative are the cumulative dose distribution data points as a function of dose (dose on x-axis).
%
% Extra parameters passed to OED will be passed onto the integrand functions.
% For example, to pass the organ specific sterilisation parameter alpha to the LinExp function,
% call OED as follows: y = OED('LinExp', dataPoints, alpha);
%
% The result is the integrated ODE dose.

%FIXME using Gaussian Quadrature method seems to be numerically unstable. Using trapezoidal method instead.
% dose=quad(@(x) LNT(x,doseCumulative),0,1);
h = 1e-5;  % <= step size
x = 0:h:1;
if strcmp(responseModel, 'LNT')
	y = LNT(x, doseCumulative);
elseif strcmp(responseModel, 'PlateauHall')
	y = PlateauHall(x, doseCumulative, varargin{:});
elseif strcmp(responseModel, 'LinExp')
	y = LinExp(x, doseCumulative, varargin{:});
else
	error('Unknown model type "%s".', responseModel);
end
dose=trapz(x, y);
return;


function y = LNT(x, doseCumulative)
y=doseInterpolate(x,doseCumulative);
return;


function y = PlateauHall(x, doseCumulative, varargin)
% Check for optional threshold parameter, otherwise set it to 4 Gy.
threshold = 4;
if length(varargin) > 0
	threshold = varargin{1};
end
d=doseInterpolate(x,doseCumulative);
y = d .* (d < threshold) + threshold .* (d >= threshold);
return;


function y = LinExp(x, doseCumulative, alpha)
d=doseInterpolate(x,doseCumulative);
y = d.*exp(-alpha.*d);
return;

