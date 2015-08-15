function result = doseIntegrate(responseModel, doseCumulative, varargin)
%result = doseIntegrate(responseModel, doseCumulative, ...)
%result = doseIntegrate(responseModel, doseCumulative, options, ...)
%
% Integrates a response model function applied to a set of dose volume histogram
% data points. i.e the following integration is performed:
%
%         1
%        /
%  F =  |  R(D(v)) dv
%       /
%      0
%
% Where R is the response model function and D is the dose as a function of
% volume fraction v. D(v) is calculated by interpolating the dose volume
% histogram data points and flipping the axes, since the data points have
% dose along the x-axis and volume along the y-axis.
%
%Parameters:
% responseModel - is the response model function that takes a dose value as the
% first parameter. Any extra parameters passed to the doseIntegrate function
% will be passed as is to this function as additional parameters. For example,
% the following call:
%   y = doseIntegrate(f, dvh, 1, 2)
% will result in the following invocation f(x, 1, 2), where x is a dose value.
%
% doseCumulative - is the cumulative dose distribution data points as a function
% of dose, with dose on the x-axis.
%
% options - is a structure with the following fields that are passed to the
% underlying integration and interpolation routines. It should be created with
% the struct() function as follows (with default values indicated):
%   options = struct('integration_method', 'quadv',
%                    'integration_tolerance', 1e-6,
%                    'interpolation_method', 'pchip');
% The structure can be empty. The option definitions are as follows:
%   'integration_method' - The integration method to use, can be one of 'quad',
%                          'quadv', 'quadl', 'quadgk' or 'quadv'.
%                          The default value used is 'trapz'.
%   'integration_tolerance' - The error tolerance parameter to use for the
%                             integration. The default tolerance is 1e-6.
%   'interpolation_method' - The interpolation method to use for the function
%                            doseInterpolate. The default method is 'pchip'.
%
% The output is the calculated integral value.
%
%Example:
% xp = 0:0.1:50;
% yp = (1+erf(10-xp))/2;
% dvh = [xp; yp];
% opts = struct('integration_method', 'trapz');
% alpha = 0.1;
% f = @(x, a) x .* exp(-a.*x);
% factor = doseIntegrate(f, dvh, opts, alpha);
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

% Authors: Artur Szostak <artursz@iafrica.com>

if nargin == 0
    % Print help message if no arguments are given.
    help doseIntegrate;
    return;
end

% Unpack the options if there are any. Otherwise the defaults are used.
integrationMethod = 'quadv';
tolerance = 1e-6;
interpMethod = 'pchip';
if length(varargin) > 0 && isstruct(varargin{1})
    integrandParams = varargin(2:length(varargin));
    opts = varargin{1};
    if isfield(opts, 'integration_method')
        integrationMethod = opts.integration_method;
    end
    if isfield(opts, 'integration_tolerance')
        tolerance = opts.integration_tolerance;
    end
    if isfield(opts, 'interpolation_method')
        interpMethod = opts.interpolation_method;
    end
else
    integrandParams = varargin;
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

integrand = @(x) invokeResponseModel(responseModel, x, doseCumulative,
                                     interpMethod, integrandParams{:});
result = integrate(integrand, 0, 1);
return;


function y = invokeResponseModel(responseModel, x, doseCumulative, interpMethod, varargin)
% This function performs an interpolation of the dose volume histogram data
% points to get dose values at certain volume fractions, before passing the
% doses to the response model function.
d = doseInterpolate(x, doseCumulative, interpMethod);
y = responseModel(d, varargin{:});
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
