function factor = RelativeRisk(model, dvh1, dvh2, varargin)
%factor = RelativeRisk(model, dvh1, dvh2, ...)
%factor = RelativeRisk(model, dvh1, dvh2, opts, ...)
%factor = RelativeRisk('LinearQuad', dvh1, dvh2, n1, n2, alpha, beta, RBEmin, RBEmax)
%factor = RelativeRisk('LinearQuad', dvh1, dvh2, opts, n1, n2, alpha, beta, RBEmin, RBEmax)
%factor = RelativeRisk('LinearQuadMultiRBE', dvh1, dvh2, n1, n2, alpha, beta, RBEmin1, RBEmax1, RBEmin2, RBEmax2)
%factor = RelativeRisk('LinearQuadMultiRBE', dvh1, dvh2, opts, n1, n2, alpha, beta, RBEmin1, RBEmax1, RBEmin2, RBEmax2)
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
% model - is the response model and can be any of: 'LinearQuad'
%
% dvh1 and dvh2 are the cumulative dose distribution data points as a function
% of dose (dose on x-axis).
%
% opts - is a structure passed to the doseIntegrate function if given. Refer to
% that function for more details.
%
% Extra parameters passed to RelativeRisk will be passed onto the integrand
% functions R1 and R2. The order and meaning depends on the value of 'model'.
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
% n2 = 11;
% alpha = 0.6;
% beta = 3;
% RBEmin = 2;
% RBEmax = 7;
% factor = RelativeRisk('LinearQuad', data1, data2, n1, n2, alpha, beta, RBEmin, RBEmax);
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

if length(varargin) > 0 && isstruct(varargin{1})
    opts = varargin{1};
    params = varargin(2:length(varargin));
else
    params = varargin;
    opts = struct();
end

% Select the integrand function based on the response model.
switch model
    case 'LinearQuad'
        n1 = params{1};
        n2 = params{2};
        alpha = params{3};
        beta = params{4};
        RBEmin = params{5};
        RBEmax = params{6};
        integrand1 = @(x) LinearQuad(x, n1, alpha, beta);
        integrand2 = @(x) LinearQuad(x, n2, RBEmax .* alpha, RBEmin.^2 .* beta);
    case 'LinearQuadMultiRBE'
        n1 = params{1};
        n2 = params{2};
        alpha = params{3};
        beta = params{4};
        RBEmin1 = params{5};
        RBEmax1 = params{6};
        RBEmin2 = params{7};
        RBEmax2 = params{8};
        integrand1 = @(x) LinearQuad(x, n1, alpha, beta);
        integrand2 = @(x) LinearQuadMulti(x, n2,
                                          RBEmax1 .* alpha, RBEmin1.^2 .* beta,
                                          RBEmax2 .* alpha, RBEmin2.^2 .* beta);
    otherwise
        error('Unknown response model type "%s".', model);
end

r1 = doseIntegrate(integrand1, dvh1, opts);
r2 = doseIntegrate(integrand2, dvh2, opts);
factor = r1 ./ r2;
return;


function y = LinearQuad(d, n, alpha, beta)
% Linear quadratic form:
%   (alpha * D + beta * D^2 / n) * exp( - (alpha * D + beta * D^2 / n) )
% where D is dose, n is the number of fractionations, alpha and beta are the
% radio-sensitivity parameters.
k = alpha .* d + beta .* d.^2 ./ n;
y = k .* exp(-k);
return;


function y = LinearQuadMulti(d, n, alpha1, beta1, alpha2, beta2)
% Linear quadratic form:
%   (alpha1 * D + beta1 * D^2 / n) * exp( - (alpha2 * D + beta2 * D^2 / n) )
% where D is dose, n is the number of fractionations, alpha and beta are the
% radio-sensitivity parameters.
k1 = alpha1 .* d + beta1 .* d.^2 ./ n;
k2 = alpha2 .* d + beta2 .* d.^2 ./ n;
y = k1 .* exp(-k2);
return;
