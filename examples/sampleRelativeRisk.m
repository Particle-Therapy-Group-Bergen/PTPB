function [results, samples] = sampleRelativeRisk(Nsamples, file1, file2, organ)
%function [results, samples] = sampleRelativeRisk(Nsamples, file1, file2, organ)
%
% Example of sampling the RelativeRisk function.
% Need to provide the following parameters:
%   Nsamples - the number of times to sample the input distributions.
%   file1 - the file from which the first (numerator) DVH will be loaded.
%   file2 - the file from which the second (denominator) DVH will be loaded.
%   organ - the organ for which to load the DVH.
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
    help sampleRelativeRisk;
    return;
end

if ~ exist('file1')
    error('Must provide a value for file1.');
end
if ~ exist('file2')
    error('Must provide a value for file2.');
end
if ~ exist('organ')
    error('Must provide a value for organ.');
end

dvh1 = getDoseVolumeHistogram(file1, organ).(organ);
dvh2 = getDoseVolumeHistogram(file2, organ).(organ);

dvh1_distrib = struct('type', 'box', 'params', {{dvh1.range_low, dvh1.range_high}});

dvh2_distrib = struct('type', 'box', 'params', {{dvh2.range_low, dvh2.range_high}});

integration_methods = struct('type', 'vector', 'params', {{'trapz'}});
integration_tolerances = struct('type', 'vector', 'params', [1e-5]);
interpolation_methods = struct('type', 'vector', 'params', {{'pchip'}});
params = struct('integration_method', integration_methods,
                'integration_tolerance', integration_tolerances,
                'interpolation_method', interpolation_methods);
opts_distrib = struct('type', 'struct', 'params', params);

% number of fractions 1.
n1_distrib = struct('type', 'vector', 'params', [1]);

% number of fractions 2.
n2_distrib = struct('type', 'vector', 'params', [1]);

% alpha radio-sensitivity parameter.
alpha_distrib = struct('type', 'gaus', 'params', {{0.1, 0.03}});

% beta radio-sensitivity parameter.
beta_distrib = struct('type', 'gaus', 'params', {{0.03, 0.005}});

% Relative Biological Effectiveness, minimum.
RBEmin_distrib = struct('type', 'triangle', 'params', {{1, 1.25, 1.5}});

% Relative Biological Effectiveness, maximum.
RBEmax_distrib = struct('type', 'triangle', 'params', {{5, 6, 7}});

[results, samples] = sampleFunction(@RR, Nsamples, dvh1_distrib, dvh2_distrib,
                                    opts_distrib, n1_distrib, n2_distrib,
                                    alpha_distrib, beta_distrib,
                                    RBEmin_distrib, RBEmax_distrib);
hist(cell2mat(results), 30);
return;


function result = RR(dvh1, dvh2, opts, n1, n2, alpha, beta, RBEmin, RBEmax)
fprintf(stdout,
        'n1 = %g\tn2 = %g\talpha = %g\tbeta = %g\tRBEmin = %g\tRBEmax = %g\n',
        n1, n2, alpha, beta, RBEmin, RBEmax);
fflush(stdout);
result = RelativeRisk('LinearQuad', dvh1, dvh2, opts, n1, n2, alpha, beta, RBEmin, RBEmax);
return;
