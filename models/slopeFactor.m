function result = slopeFactor(cancerSite, gender, exposureAge, age)
%result = slopeFactor(cancerSite, gender, exposureAge, age)
%
% Calculates the slope factor that must be multiplied with the average dose or organ equivalent dose to estimate the risk.
%
%Where,
% cancerSite selects the parameters for a specified organ.
%
% gender selects the gender of the patient.
%
% exposureAge is the age of patient at exposure.
%
% age is the attained age.
%
% The output is the coefficient that must be multiplied with dose to get the secondary cancer risk.
%
%Example: a = 10:40; plot(a, slopeFactor('BEIR_EAR_Lung', 'F', 10, a));
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
    help slopeFactor;
    return;
end

[beta, gamma, eta, tau, c1, c2, c3, c4, k1, k2] = selectParams(cancerSite, gender);
% Calculating:
%   beta * exp(gamma * e_star + eta * ln(a_star))
%
% where
%   e_star = { c1 * e + c2  if e < tau
%            { c3 * e + c4  if e >= tau
%
%   a_star = k1 * a + k2
%
%   e = age at exposure
%   a = attained age
%   tau = threshold parameter
%   c1..4 and k1..2 are affine transformation parameters.
%
% tau, c1..4, k1..2 etc are chosen to match appropriate parametric fits from different reports, eg. BEIR VII.
e_star = (c1 .* exposureAge + c2).*(exposureAge < tau) + (c3 .* exposureAge + c4).*(exposureAge >= tau);
a_star = k1 .* age + k2;
result = beta .* exp(gamma .* e_star + eta .* log(a_star));
return;


function [beta, gamma, eta, tau, c1, c2, c3, c4, k1, k2] = selectParams(cancerSite, gender)

% Update data table as needed:
% BEIR VII: Table 12-2 Committee's preferred ERR and EAR models for estimating site-specific solid cancer incidence and
% mortality. Some parameters estimated by Gonzales et al using RadRAT.
%                                            beta_M  beta_F  gamma   eta  tau  c1    c2  c3  c4  k1    k2
table = struct('BEIR_ERR_Stomach',          [0.21    0.48    -0.3   -1.4  30   1/10  -3  0   0   1/60  0  ],
               'BEIR_ERR_Colon',            [0.63    0.43    -0.3   -1.4  30   1/10  -3  0   0   1/60  0  ],
               'BEIR_ERR_Liver',            [0.32    0.32    -0.3   -1.4  30   1/10  -3  0   0   1/60  0  ],
               'BEIR_ERR_Lung',             [0.32    1.4     -0.3   -1.4  30   1/10  -3  0   0   1/60  0  ],
               'BEIR_ERR_Breast',           [0       0.51     0     -2     0   1      0  1   0   1/60  0  ],
               'BEIR_ERR_Prostate',         [0.12    0       -0.3   -1.4  30   1/10  -3  0   0   1/60  0  ],
               'BEIR_ERR_Uterus',           [0       0.055   -0.3   -1.4  30   1/10  -3  0   0   1/60  0  ],
               'BEIR_ERR_Ovary',            [0       0.38    -0.3   -1.4  30   1/10  -3  0   0   1/60  0  ],
               'BEIR_ERR_Bladder',          [0.50    1.65    -0.3   -1.4  30   1/10  -3  0   0   1/60  0  ],
               'BEIR_ERR_Thyroid',          [0.53    1.05    -0.83   0     0   1      0  1   0   1/60  0  ],
               'RadRAT_ERR_Rectum',         [0.12    0.12    -0.3   -1.4  30   1/10  -3  0   0   1/60  0  ],
               'BEIR_ERR_OtherSolidCancer', [0.27    0.45    -0.3   -2.8  30   1/10  -3  0   0   1/60  0  ],

               'BEIR_EAR_Stomach',          [4.9     4.9     -0.41   2.8  30   1/10  -3  0   0   1/60  0  ],
               'BEIR_EAR_Colon',            [3.2     1.6     -0.41   2.8  30   1/10  -3  0   0   1/60  0  ],
               'BEIR_EAR_Liver',            [2.2     1.0     -0.41   4.1  30   1/10  -3  0   0   1/60  0  ],
               'BEIR_EAR_Lung',             [2.3     3.4     -0.41   5.2  30   1/10  -3  0   0   1/60  0  ],
               'BEIR_EAR_Breast',           [0       9.4     -0.51   3.5   0   1      0  1   0   1/60  0  ],
               'BEIR_EAR_Prostate',         [0.11    0       -0.41   2.8  30   1/10  -3  0   0   1/60  0  ],
               'BEIR_EAR_Uterus',           [0       1.2     -0.41   2.8  30   1/10  -3  0   0   1/60  0  ],
               'BEIR_EAR_Ovary',            [0       0.7     -0.41   2.8  30   1/10  -3  0   0   1/60  0  ],
               'BEIR_EAR_Bladder',          [1.2     0.75    -0.41   6.0  30   1/10  -3  0   0   1/60  0  ],
               'BEIR_EAR_Thyroid',          [0       0        0      0     0   0      0  0   0    0    0  ],
               'RadRAT_Rectum',             [0.34    0.34    -0.41   2.8  30   1/10  -3  0   0   1/60  0  ],
               'BEIR_EAR_OtherSolidCancer', [6.2     4.8     -0.41   2.8  30   1/10  -3  0   0   1/60  0  ]); 

if ~ isfield(table, cancerSite)
    error('Invalid cancerSite value "%s".', cancerSite);
end
row = getfield(table, cancerSite);
switch gender
    case 'M'
        beta = row(1);
    case 'F'
        beta = row(2);
    otherwise
        error('Invalid gender string. Must be one of "F" or "M".');
end
gamma = row(3);
eta = row(4);
tau = row(5);
c1 = row(6);
c2 = row(7);
c3 = row(8);
c4 = row(9);
k1 = row(10);
k2 = row(11);
return;
