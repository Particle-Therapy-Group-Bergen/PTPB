function factor = RelativeRisk(options, responseModel, doseCumulative1, doseCumulative2, varargin)
%factor = RelativeRisk(options, responseModel, doseCumulative1, doseCumulative2, ...)
%
% Calculates the Relative Risk (RR) from the cumulative dose distribution data points.
% i.e the following integration is performed:
%
%         1,1                                  1,1
%         /                                    /    R1( D1(v1) )
%  RR =   |  R( D1(v1), D2(v2) )  dv1 dv2  =   |    ------------  dv1 dv2
%        /                                    /     R2( D2(v2) )
%       0,0                                  0,0
%
% Where R is the relative risk function, which itself is a ratio of two functions
% based on come response models. R1 and R2 are the probability of malignant
% transformations per cell.
% D1 and D2 are each the dose as a function of volume fraction corresponding to
% R1 and R2 respectively. Di(vi) is calculated by interpolating the cumulative
% dose distribution data points and flipping the axes, since the data points have
% dose along the x-axis and volume along the y-axis.
%
%Parameters:
% options is a structure with the following parameters that are passed to the
% underlying integration and interpolation routines. It should be created with
% the struct() function as follows (with default values indicated):
%   options = struct('integration_method', 'quadv',
%                    'integration_tolerance', 1e-6,
%                    'interpolation_method', 'pchip');
% The structure can be empty. The option definitions are as follows:
%   'integration_method' - The integration method to use, can be one of 'quad',
%                          'quadv', 'quadl', 'quadgk' or 'trapz'.
%   'integration_tolerance' - The error tolerance parameter to use for the integration.
%   'interpolation_method' - The interpolation method to use for doseInterpolate.
%
% responseModel can be any of: 'LinearQuad'
%
% doseCumulative1 and doseCumulative1 are the cumulative dose distribution data
% points as a function of dose (dose on x-axis).
%
% Extra parameters passed to RelativeRisk will be passed onto the integrand functions.
% For example, to pass the organ specific alpha and beta parameters for the
% 'LinearQuad' method, call RelativeRisk as follows:
%   y = RelativeRisk('LinExp', dataPoints, alpha);
%
% The output is the integrated relative risk factor.
%
%Example:
% xp = 0:0.1:50;
% yp = (1+erf(10-xp))/2;
% data1 = [xp; yp];
% xp = 0:0.1:50;
% yp = (1+erf(7-xp*0.5))/2;
% data2 = [xp; yp];
% n1 = 5;
% alpha1 = 0.5;
% beta1 = 2;
% n2 = 11;
% alpha2 = 0.6;
% beta2 = 3;
% RBEmin = 2;
% RBEmax = 7;
% factor = RelativeRisk(struct('integration_method', 'trapz'), 'LinearQuad',
%                       data1, data2, n1, alpha1, beta1, n2, alpha2, beta2,
%                       RBEmin, RBEmax);
%

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

% Authors: Artur Szostak <artursz@iafrica.com>, Camilla H Stokkevaag <camilla.stokkevag@ift.uib.no>

if nargin == 0
    % Print help message if no arguments are given.
    help RelativeRisk;
    return;
end

% Unpack the options if there are any. Otherwise the defaults are used.
integrationMethod = 'trapz';
tolerance = 1e-6;
interpMethod = 'pchip';
if length(options) ~= 0
    if isstruct(options)
        if isfield(options, 'integration_method')
            integrationMethod = options.integration_method;
        end
        if isfield(options, 'integration_tolerance')
            tolerance = options.integration_tolerance;
        end
        if isfield(options, 'interpolation_method')
            interpMethod = options.interpolation_method;
        end
    else
        error('The options field must be a structure.');
    end
end

% Select the integrand function based on the response model.
switch responseModel
    case 'LinearQuad'
        n1 = varargin{1};
        alpha1 = varargin{2};
        beta1 = varargin{3};
        n2 = varargin{4};
        alpha2 = varargin{5};
        beta2 = varargin{6};
        RBEmin = varargin{7};
        RBEmax = varargin{8};
        integrand1 = @(x) LinearQuad(x, interpMethod, doseCumulative1,
                                     n1, alpha1, beta1);
        integrand2 = @(x) LinearQuad(x, interpMethod, doseCumulative2,
                                     n2, RBEmax .* alpha2, RBEmin.^2 .* beta2);
    otherwise
        error('Unknown response model type "%s".', responseModel);
end

% Select the integration function based on the integration method.
switch integrationMethod
    case 'quad'
        integrate = @(f, a, b) quad(f, a, b, tolerance);
    case 'quadv'
        integrate = @(f, a, b) quadv(f, a, b, tolerance);
    case 'quadl'
        integrate = @(f, a, b) quadl(f, a, b, tolerance);
    case 'quadgk'
        integrate = @(f, a, b) quadgk(f, a, b, tolerance);
    case 'trapz'
        integrate = @(f, a, b) trapzIntegrate(f, a, b, tolerance);
    otherwise
        error('Unsupported integration method "%s".', integrationMethod);
end

factor = integrate(integrand1, 0, 1) / integrate(integrand2, 0, 1);
return;


function y = trapzIntegrate(f, a, b, tol)
% Simple integration method using the trapz function.

h = abs(b-a)*tol*2;  % <= initial step size
x = a:h:b;
oldy = trapz(x, f(x));
% Halve the step size and calculate again:
h = h/2;
x = a:h:b;
y = trapz(x, f(x));
% While the difference (estimate of error) is greater than the tolerance
% threshold, keep halving the step size and calculate again.
while abs(y - oldy) > tol
    oldy = y;
    h = h/2;
    x = a:h:b;
    y = trapz(x, f(x));
end
return;


function y = LinearQuad(x, interpMethod, doseCumulative, n, alpha, beta)
% Linear quadratic form:
%   (alpha * D + beta * D^2 / n) * exp( - (alpha * D + beta * D^2 / n) )
d = doseInterpolate(x, doseCumulative, interpMethod);
k = alpha .* d + beta .* d.^2 ./ n;
y = k .* exp(-k);
return;
