function makeToleranceScanPlot(X, varargin)
% Makes a plot of integration performance as a function of tolerance.

Y = [];

for n = 1:length(X)
    data = load(varargin{n}, 'Results');
    Y(:,n) = data.Results(:,1);
    if ~ exist('alpha')
        alpha = data.Results(1, 2);
        beta = data.Results(1, 3);
        RBEmin = data.Results(1, 4);
        RBEmax = data.Results(1, 5);
    end
end

semilogxerr(X, mean(Y), std(Y), 'x');
%loglogerr(X, mean(Y), std(Y), 'x');
title(sprintf('alpha = %g, beta = %g, RBEmin = %g, RBEmax = %g', alpha, beta, RBEmin, RBEmax));
xlabel('Integration tolerance');
ylabel('Relative Risk');
