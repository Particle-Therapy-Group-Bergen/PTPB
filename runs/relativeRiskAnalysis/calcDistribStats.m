function [results, uncertainty, datamat] = calcDistribStats(filename, M, N)
%[results, uncertainty, datamat] = calcDistribStats(filename, M, N)
%
%Calculates various statistical parameters for the sample data loaded from the
%given file.
%
%filename - The name of the file to load sample data from.
%M - The number of boot-strap samples to create to estimate the uncertainty of
%    the statistical parameters.
%N - The number of input samples to compute the parameters from.
%
%The following outputs are produced:
%
%results - The values of the statistical parameters. The vector of values
%          include: [Ql, median, Qh, mean, std.dev., Qh-Ql]
%          Where Ql is the lower 2.5% quantile and Qh is the upper 97.5%
%          quantile.
%uncertainty - The uncertainties of the statistical parameters.
%datamat - The statistical parameters samples used to estimate the
%          uncertainties of the final results.

if ~ exist('filename')
    filename = 'cion_bladder_sample_data.mat';
end
if ~ exist('M')
    M = 100;
end

load(filename);
if ~ exist('N')
    N = size(Results, 1);
end

datamat = zeros(M,5);
for m = 1:M
    %D = Results(1:N,1);   % not removing outliers
    D = Results( find(0 < Results(1:N,1) & Results(1:N,1) < 100) , 1);   % removing outliers
    S = bootStrap(D, 1, 'random')';
    Q = quantile(S, [0.025, 0.5, 0.975], 1, 8)';
    datamat(m,1:3) = Q;
    datamat(m,4) = mean(S);
    datamat(m,5) = std(S);
    datamat(m,6) = Q(3) - Q(1);
end
results(1,:) = mean(datamat);
uncertainty(1,:) = std(datamat);
