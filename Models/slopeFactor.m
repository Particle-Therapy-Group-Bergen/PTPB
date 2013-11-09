function result = slopeFactor(cancerSite, gender, exposureAge, age)
% Calculates the slope factor that must be multiplied with the average dose or organ equivalent dose to estimate the risk.
% cancerSite selects the parameters for a specified organ.
% gender selects the gender of the patient.
% exposureAge is the age of patient at exposure.
% age is the attained age.

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
%                                  beta_M  beta_F  gamma  eta   tau  c1    c2  c3  c4  k1    k2
table = struct('BEIR_VII_Lung',   [0.32    1.4     -0.3   -1.4  30   1/10  -3  0   0   1/60  0  ],
               'BEIR_VII_Breast', [0       0.51     0     -2    0    1      0  1   0   1/60  0  ]);

if ! isfield(table, cancerSite)
	error('Invalid cancerSite value "%s".', cancerSite);
end
row = getfield(table, cancerSite);
if strcmp(gender, 'M')
	beta = row(1);
elseif strcmp(gender, 'F')
	beta = row(2);
else
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

