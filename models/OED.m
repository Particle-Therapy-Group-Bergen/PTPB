function dose = OED(model, dvh, varargin)
%dose = OED(model, dvh, ...)
%dose = OED(model, dvh, opts, ...)
%dose = OED('LNT', dvh)
%dose = OED('LNT', dvh, opts)
%dose = OED('PlateauHall', dvh, threshold)
%dose = OED('PlateauHall', dvh, opts, threshold)
%dose = OED('LinExp', dvh, alpha)
%dose = OED('LinExp', dvh, opts, alpha)
%dose = OED('Competition', dvh, alpha1, beta1, alpha2, beta2, n)
%dose = OED('Competition', dvh, opts, alpha1, beta1, alpha2, beta2, n)
%dose = OED('LinPlat', dvh, delta)
%dose = OED('LinPlat', dvh, opts, delta)
%
% Calculates the Organ Equivalent Dose (OED) from the cumulative dose volume
% histogram data points. i.e the following integration is performed:
%
%          1
%         /
%  OED =  | R( D(v) ) dv
%        /
%        0
%
% Where R is the dose response model and D is the dose as a function of volume
% fraction. D(v) is calculated by interpolating the cumulative dose volume
% histogram data points and flipping the axes, since the data points have dose
% along the x-axis and volume along the y-axis.
%
%Parameters:
% model - is the response model name and can be any of: 'LNT', 'PlateauHall',
% 'LinExp', 'Competition' or 'LinPlat'.
%
% dvh - are the cumulative dose volume histogram data points as a function of
% dose (dose on x-axis).
%
% opts - is a structure passed to the doseIntegrate function if given. Refer to
% that function for more details.
%
% Extra parameters passed to OED will be passed onto the integrand functions.
% For example, to pass the organ specific sterilisation parameter alpha to the
% LinExp function, call OED as follows: y = OED('LinExp', dataPoints, alpha);
%
% The output is the integrated OED dose.
%
%Example:
% xp = 0:0.1:50;
% yp = (1+erf(10-xp))/2;
% data = [xp; yp];
% threshold = 35;
% dose = OED('PlateauHall', data, struct('integration_method', 'trapz'), threshold);
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
    help OED;
    return;
end

% Select the integrand function based on the response model and integrate.
switch model
    case 'LNT'
        dose = doseIntegrate(@LNT, dvh, varargin{:});
    case 'PlateauHall'
        dose = doseIntegrate(@PlateauHall, dvh, varargin{:});
    case 'LinExp'
        dose = doseIntegrate(@LinExp, dvh, varargin{:});
    case 'Competition'
        dose = doseIntegrate(@Competition, dvh, varargin{:});
    case 'LinPlat'
        dose = doseIntegrate(@LinPlat, dvh, varargin{:});
    otherwise
        error('Unknown response model type "%s".', model);
end
return;


function y = LNT(d)
y = d;
return;


function y = PlateauHall(d, threshold)
% Check for optional threshold parameter, otherwise set it to 4 Gy.
if ~ exist('threshold')
    threshold = 4.5;
end
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
