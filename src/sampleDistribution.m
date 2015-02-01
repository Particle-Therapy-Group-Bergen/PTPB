function result = sampleDistribution(distribution, N, varargin)
%result = sampleDistribution(distribution, N, ...)
%result = sampleDistribution('delta', N, value)
%result = sampleDistribution('double_delta', N, A, B)
%result = sampleDistribution('box', N, min, max)
%result = sampleDistribution('box95', N, lower, upper)
%result = sampleDistribution('triangle', N, min, mode, max)
%result = sampleDistribution('triangle95', N, lower, mean, upper)
%result = sampleDistribution('gaus', N, mu, sigma)
%result = sampleDistribution('gaus95', N, lower, upper)
%result = sampleDistribution('lognorm', N, mu, sigma)
%result = sampleDistribution('lognorm95', N, lower, upper)
%
% Randomly samples a distribution.
%
%Where,
% distribution is a string naming the distribution type to sample. This can be
% one of the following values:
%   'delta' - A delta distribution requires one parameter indicating the value
%       to returun.
%   'double_delta' - A distribution with two delta values (two states), requires
%       two parameters A and B, one for each value.
%   'box' - The box distribution requires two parameters, the lower bound of
%       the distribution and the upper bound of the distribution.
%   'box95' - A box distribution defined by a pair of 95% confidence interval
%       (CI) parameters, the lower 95% CI bound and upper CI bound.
%   'triangle' - The triangle distribution requires 3 parameters, the lower
%       bound of the distribution, the mode (peak) of the distribution and the
%       upper bound B.
%   'triangle95' - A triangle distribution defined by 3 parameters, the lower
%       bound of the 95% CI interval, the mean of the distribution and the upper
%       bound of the CI interval.
%   'triangle95mode' - A triangle distribution similar defined by 3 parameters,
%       the lower bound of the 95% CI interval, the mode (peak) of the
%       distribution and the upper bound of the CI interval.
%   'gaus' - A gaussian distribution requires two parameters, the mean (mu) and
%       standard deviation (sigma).
%   'gaus95' - A gaussian distribution defined by two parameters, the lower and
%       upper bound of the CI interval.
%   'lognorm' - The log normal distribution requires two parameters, the mean
%       and sigma of the log of the distribution.
%   'lognorm95' - A log normal distribution defined by requires two parameters,
%       the lower and upper bound of the CI interval.
%
% N is the number of samples to produce.
%
% All additional parameters are used to define the shape of the distribution
% being sampled. These can be vectors, as long as they are the same length.

if nargin == 0
    % Print help message if no arguments are given.
    help sampleDistribution;
    return;
end
if ~ exist('N')
    error('Number of samples to generate N was not given.')
end

switch distribution
    case 'delta'
        if nargin < 3
            error('Require at least one distribution parameter.')
        end
        result = repmat(varargin{1}, N, 1);
    case 'double_delta'
        if nargin < 4
            error('Require at least two distribution parameters.')
        end
        A = repmat(varargin{1}, N, 1);
        B = repmat(varargin{2}, N, 1);
        result = rand2state(A, B);
    case 'box'
        if nargin < 4
            error('Require at least two distribution parameters.')
        end
        a = min(varargin{2}, varargin{1});
        b = max(varargin{2}, varargin{1});
        A = repmat(a, N, 1);
        B = repmat(b, N, 1);
        result = randBox(A, B);
    case 'box95'
        if nargin < 4
            error('Require at least two distribution parameters.')
        end
        a = min(varargin{2}, varargin{1});
        b = max(varargin{2}, varargin{1});
        delta = (b - a) * (1 / 0.95 - 1) / 2;
        A = repmat(a - delta, N, 1);
        B = repmat(b + delta, N, 1);
        result = randBox(A, B);
    case 'triangle'
        if nargin < 5
            error('Require at least three distribution parameters.')
        end
        A = repmat(varargin{1}, N, 1);
        B = repmat(varargin{3}, N, 1);
        C = repmat(varargin{2}, N, 1);
        result = randTriangle(A, B, C);
    case 'triangle95'
        if nargin < 5
            error('Require at least three distribution parameters.')
        end
        l = varargin{1};
        mu = varargin{2};
        u = varargin{3};
        Cl = 0.025;
        Cu = 0.975;
        % Find an initial approximation for a and b.
        c = 3*mu-l-u;
        g = @(x) triangleEquationsMode(x, l, c, u, Cl, Cu);
        a = c - (c-l)/(1-Cl);
        b = (u-c)/Cu + c;
        r = fsolve(g, [a, b], struct('TolX', 1e-12, 'TolFun', 1e-12));
        a = r(1);
        b = r(2);
        % Now try find exact values for a, b, and c.
        f = @(x) triangleEquations(x, l, mu, u, Cl, Cu);
        r = fsolve(f, [a, b, c], struct('TolX', 1e-12, 'TolFun', 1e-12));
        a = r(1);
        b = r(2);
        c = r(3);
        A = repmat(a, N, 1);
        B = repmat(b, N, 1);
        C = repmat(c, N, 1);
        result = randTriangle(A, B, C);
    case 'triangle95mode'
        if nargin < 5
            error('Require at least three distribution parameters.')
        end
        l = varargin{1};
        c = varargin{2};
        u = varargin{3};
        Cl = 0.025;
        Cu = 0.975;
        f = @(x) triangleEquationsMode(x, l, c, u, Cl, Cu);
        a = c - (c-l)/(1-Cl);
        b = (u-c)/Cu+c;
        r = fsolve(f, [a, b], struct('TolX', 1e-12, 'TolFun', 1e-12));
        a = r(1);
        b = r(2);
        A = repmat(a, N, 1);
        B = repmat(b, N, 1);
        C = repmat(c, N, 1);
        result = randTriangle(A, B, C);
    case 'gaus'
        if nargin < 4
            error('Require at least two distribution parameters.')
        end
        mu = repmat(varargin{1}, N, 1);
        sigma = repmat(varargin{2}, N, 1);
        result = randn(size(mu)) .* sigma + mu;
    case 'gaus95'
        if nargin < 4
            error('Require at least two distribution parameters.')
        end
        l = repmat(varargin{1}, N, 1);
        u = repmat(varargin{2}, N, 1);
        Cl = 0.025;
        Cu = 0.975;
        k1 = sqrt(2).*erfinv(2.*Cl - 1);
        k2 = sqrt(2).*erfinv(2.*Cu - 1);
        mu = (k2.*l - k1.*u) ./ (k2 - k1);
        sigma = (l-mu) / k1;
        result = randn(size(mu)) .* sigma + mu;
    case 'lognorm'
        if nargin < 4
            error('Require at least two distribution parameters.')
        end
        mu = repmat(varargin{1}, N, 1);
        sigma = repmat(varargin{2}, N, 1);
        result = exp(randn(size(mu)) .* sigma + mu);
    case 'lognorm95'
        if nargin < 4
            error('Require at least two distribution parameters.')
        end
        l = repmat(varargin{1}, N, 1);
        u = repmat(varargin{2}, N, 1);
        Cl = 0.025;
        Cu = 0.975;
        k1 = sqrt(2).*erfinv(2.*Cl - 1);
        k2 = sqrt(2).*erfinv(2.*Cu - 1);
        mu = (k2.*log(l) - k1.*log(u)) ./ (k2 - k1);
        sigma = (log(l)-mu) / k1;
        result = exp(randn(size(mu)) .* sigma + mu);
    otherwise
        error('Distribution type "%s" is not supported.', distribution);
end
return;


function X = randBox(a, b)
% Generate random number matrix with the same dimentions as a, distributed according to
% a box distribution with min a, max b.

U = rand(size(a));
X = (b-a).*U + a;
return;


function X = randTriangle(a, b, c)
% Generate random number matrix with the same dimentions as a, distributed according to
% a triangular distribution with min a, max b and mode (peak) c.

U = rand(size(a));
F = (b-a).*U < (c-a);
X = (a + sqrt(U.*(b-a).*(c-a))) .* F + (b - sqrt((1-U).*(b-a).*(b-c))).*(~F);
return;


function X = rand2state(a, b)
% Generate random number matrix where each entry can have only one of 2 values: a or b.

U = round( rand(size(a)) );
X = a.*U + b.*(1-U);
return;


function Y = triangleEquations(X, l, mu, u, Cl, Cu)
% System of linear equations to solve to find a, b, and c parameters given the
% confidence interval [l .. u], the mean mu, and the quantile values [Cl, Cu]
% that define the 95% confidence interval.
%  3 * mu - a - b - a = 0
%  (l - a)^2 / ((b - a) * (c - a)) - Cl = 0
%  (b - u)^2 / ((b - a) * (b - c)) + Cu - 1 = 0

a = X(1);
b = X(2);
c = X(3);
Y(1) = 3*mu - a - b - c;
Y(2) = (l-a)^2 / ((b-a)*(c-a)) - Cl;
Y(3) = (b-u)^2 / ((b-a)*(b-c)) + Cu - 1;
return;


function Y = triangleEquationsMode(X, l, c, u, Cl, Cu)
% System of linear equations to solve to find a and b parameters given the
% confidence interval [l .. u], the mode (peak) c, and the quantile values
% [Cl, Cu] that define the 95% confidence interval.
%  (l - a)^2 / ((b - a) * (c - a)) - Cl = 0
%  (b - u)^2 / ((b - a) * (b - c)) + Cu - 1 = 0

a = X(1);
b = X(2);
Y(1) = (l-a)^2 / ((b-a)*(c-a)) - Cl;
Y(2) = (b-u)^2 / ((b-a)*(b-c)) + Cu - 1;
return;
