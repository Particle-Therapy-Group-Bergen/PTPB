function relativeRiskPlots()
%function relativeRiskPlots()
%
% This example script reproduces plots from [1] figure 5. This is an example of
% using the RelativeRisk function.
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

alpha = 0.1;  % radio-sensitivity parameter.
beta = 0.03;  % radio-sensitivity parameter.
RBEmin = 1.03;  % Relative Biological Effectiveness, minimum.
RBEmax = 1.25;  % Relative Biological Effectiveness, maximum.
ProtonRBE = 1.1; % Relative Biological Effectiveness for protons.

x = 0.001:0.2:10;

figure;
y1 = RR(x, 1, alpha, beta, RBEmin, RBEmax, ProtonRBE);
y2 = RR(x, 11, alpha, beta, RBEmin, RBEmax, ProtonRBE);
y3 = RR(x, 21, alpha, beta, RBEmin, RBEmax, ProtonRBE);
y4 = RR(x, 31, alpha, beta, RBEmin, RBEmax, ProtonRBE);
y5 = RR(x, 41, alpha, beta, RBEmin, RBEmax, ProtonRBE);
plot(x, y1, x, y2, x, y3, x, y4, x, y5);
title('CoGyEq protons c.f. megavoltage x-rays');
xlabel('Dose per fraction (Gy)');
ylabel('Relative risk per cell');
legend('1 fraction', '11 fractions', '21 fractions', '31 fractions', '41 fractions');

RBEmin = 1.25;
RBEmax = 6;
CarbonRBE = 3; % Relative Biological Effectiveness for carbon ions.

figure;
y1 = RR(x, 1, alpha, beta, RBEmin, RBEmax, CarbonRBE);
y2 = RR(x, 6, alpha, beta, RBEmin, RBEmax, CarbonRBE);
y3 = RR(x, 11, alpha, beta, RBEmin, RBEmax, CarbonRBE);
y4 = RR(x, 16, alpha, beta, RBEmin, RBEmax, CarbonRBE);
y5 = RR(x, 21, alpha, beta, RBEmin, RBEmax, CarbonRBE);
plot(x, y1, x, y2, x, y3, x, y4, x, y5);
title('C Ions c.f. megavoltage x-rays');
xlabel('Dose per fraction (Gy)');
ylabel('Relative risk per cell');
legend('1 fraction', '6 fractions', '11 fractions', '16 fractions', '21 fractions');
return;


function risk = RR(dose, N, alpha, beta, RBEmin, RBEmax, scalingRBE)

% Prepare a function to produce simulated DVH data as a step function.
simdvh = @(d) [0, 1; d, 1; d, 0; 2*d, 0];

opts = struct('integration_method', 'trapz', 'integration_tolerance', 1e-3,
              'interpolation_method', 'linear');

risk = [];
for n = 1:length(dose)
    d = dose(n);
    dvh1 = simdvh(N*d);
    dvh2 = simdvh(N*d/scalingRBE);
    risk(n) = RelativeRisk('LinearQuad', dvh1, dvh2, opts, N, N, alpha, beta,
                           RBEmin, RBEmax);
end
