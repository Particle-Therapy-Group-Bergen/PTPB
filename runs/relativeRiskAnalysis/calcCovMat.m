function C = calcCovMat(filename)
%C = calcCovMat(filename)
%
% Calculates the covariance matrix for the sample data loaded from 'filename'.

load(filename);
X = Results;
n = size(X, 1);
M = repmat(mean(X), n, 1);
A = X - M;
C = A'*A / (n-1);
