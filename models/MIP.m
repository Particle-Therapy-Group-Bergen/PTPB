function result = MIP(model, dvh, varargin)
%result = MIP(model, dvh, ...)
%result = MIP(model, dvh, opts, ...)
%result = MIP('LinearQuad', dvh, n, alpha, beta)
%result = MIP('LinearQuad', dvh, opts, n, alpha, beta)
%
% Calculates the Malignancy Induction Probability (MIP) per cell from the
% cumulative dose volume histogram data points. i.e the following integration
% is performed:
%
%          1
%         /
%  MIP =  | R( D(v) ) dv
%        /
%        0
%
% Where R is the dose response model and D is the dose as a function of volume
% fraction. D(v) is calculated by interpolating the cumulative dose volume
% histogram data points and flipping the axes, since the data points have dose
% along the x-axis and volume along the y-axis.
%
%Parameters:
% model - is the response model and can be any of: 'LinearQuad'
%
% dvh - is the cumulative dose volume histogram data points as a function of
% dose, with dose on the x-axis.
%
% opts - is a structure passed to the doseIntegrate function if given. Refer to
% that function for more details.
%
% Extra parameters passed to MIP will be passed onto the response model integrand
% function.
%
% The output is the integrated MIP probability.
%
%Example:
% xp = 0:0.1:50;
% yp = (1+erf(10-xp))/2;
% data = [xp; yp];
% dose = MIP('LinearQuad', data, 5, 0.1, 0.03);
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
    help MIP;
    return;
end

% Select the integrand function based on the response model and integrate.
switch model
    case 'LinearQuad'
        dose = doseIntegrate(@LinearQuad, dvh, varargin{:});
    otherwise
        error('Unknown response model type "%s".', model);
end
return;

function y = LinearQuad(d, n, alpha, beta)
% Linear quadratic form:
%   (alpha * D + beta * D^2 / n) * exp( - (alpha * D + beta * D^2 / n) )
k = alpha .* d + beta .* d.^2 ./ n;
y = k .* exp(-k);
return;
