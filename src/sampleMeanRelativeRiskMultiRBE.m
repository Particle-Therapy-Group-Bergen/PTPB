function result = sampleMeanRelativeRiskMultiRBE(Nsamples, filepat1, filepat2,
    patients, organ, n1, n2, scale1, scale2, alpha_distrib, beta_distrib,
    RBEmin1_distrib, RBEmax1_distrib, RBEmin2_distrib, RBEmax2_distrib,
    opts, namemap)
%function result = sampleMeanRelativeRiskMulti(Nsamples, filepat1, filepat2,
%  patients, organ, n1, n2, scale1, scale2, alpha_distrib, beta_distrib,
%  RBEmin1_distrib, RBEmax1_distrib, RBEmin2_distrib, RBEmax2_distrib,
%  opts, namemap)
%
% Calculates sample points from the mean relative risk distribution, given
% patient input files, fractionation, scaling factors, alpha, beta and RBE
% parameter distributions to use.
%
% The results will be a matrix with each sample on a row and the following
% columns: relative risk result, alpha, beta, RBEmin and RBEmax.
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

if ~ exist('Nsamples')
    Nsamples = 10;
end

% Input DVH data file pattern 1.
if ~ exist('filepat1')
    filepat1 = 'data/VMATdvh/vmat%d.mat';
end

% number of fractions 1.
if ~ exist('n1')
    n1 = 25;
end

% Dose scaling factor for DVH 1. Dose will be multiplied by this.
if ~ exist('scale1')
    scale1 = 1;
end

% Input DVH data file pattern 2.
if ~ exist('filepat2')
    filepat2 = 'data/CionDataPhysicalDose/HUH%dphysical_dvh.mat';
end

% number of fractions 2.
if ~ exist('n2')
    n2 = 12;
end

% Dose scaling factor for DVH 2. Dose will be multiplied by this.
if ~ exist('scale2')
    scale2 = 1;
end

if ~ exist('patients')
    patients = [12, 33, 35, 36, 37, 39, 41, 42, 43, 44];
end

if ~ exist('organ')
    organ = 'Bladder';
end

if ~ exist('alpha_distrib')
    alpha_distrib = struct('type', 'gaus', 'params', {{0.1, 0.03}});
end
if ~ exist('beta_distrib')
    beta_distrib = struct('type', 'gaus', 'params', {{0.03, 0.005}});
end
if ~ exist('RBEmin1_distrib')
    RBEmin1_distrib = struct('type', 'triangle', 'params', {{1.1, 1.25, 1.4}});
end
if ~ exist('RBEmax1_distrib')
    RBEmax1_distrib = struct('type', 'triangle', 'params', {{5, 6, 7}});
end
if ~ exist('RBEmin2_distrib')
    RBEmin2_distrib = struct('type', 'triangle', 'params', {{1.1, 1.25, 1.4}});
end
if ~ exist('RBEmax2_distrib')
    RBEmax2_distrib = struct('type', 'triangle', 'params', {{5, 6, 7}});
end

% Integration and interpolation options.
% Note, we want 1 randomised bootstrap sampling for each single input DVH
% ampling by default. Thus, bootstrap_samples = 1:
if ~ exist('opts')
    opts = struct('integration_method', 'trapz',
                  'integration_tolerance', 1e-5,
                  'interpolation_method', 'linear',
                  'sample_dvh', 1,
                  'bootstrap_samples', 1,
                  'bootstrap_method', 'random');
end

% Organ name remapping:
namemap = {
        'Bladder_P', 'Bladder';
        'Rectum_P_MT', 'Rectum';
    };

% Load DVH data from files.
dvh1list = {};
dvh2list = {};
for k = 1:length(patients)
    % Relative Risk = file1/file2
    file1 = sprintf(filepat1, patients(k));
    file2 = sprintf(filepat2, patients(k));
    dvh1 = getDoseVolumeHistogram(file1, namemap, organ).(organ).datapoints;
    dvh1(:,1) = dvh1(:,1) .* scale1;
    dvh1list{k} = dvh1;
    dvh2 = getDoseVolumeHistogram(file2, namemap, organ).(organ).datapoints;
    dvh2(:,1) = dvh2(:,1) .* scale2;
    dvh2list{k} = dvh2;
end

if isfield(opts, 'bootstrap_samples')
    bootstrap_samples = opts.bootstrap_samples;
else
    bootstrap_samples = 1;
end
if isfield(opts, 'bootstrap_method')
    bootstrap_method = opts.bootstrap_method;
else
    bootstrap_method = 'random';
end

if isfield(opts, 'sample_dvh') && opts.sample_dvh
    % Prepare the sampling ranges for the DVHs.
    params1 = {};
    params2 = {};
    for k = 1:length(patients)
        [low, high] = getDVHBinRanges(dvh1list{k});
        params1{k} = struct('type', 'box', 'params', {{low, high}});
        [low, high] = getDVHBinRanges(dvh2list{k});
        params2{k} = struct('type', 'box', 'params', {{low, high}});
    end
    dvh1_distrib = struct('type', 'array', 'params', {params1});
    dvh2_distrib = struct('type', 'array', 'params', {params2});

    func = @(dvh1list, dvh2list, alpha, beta, RBEmin1, RBEmax1, RBEmin2, RBEmax2) CalculateMeanRR(
                dvh1list, dvh2list, opts, n1, n2, alpha, beta, RBEmin1, RBEmax1,
                RBEmin2, RBEmax2, bootstrap_samples, bootstrap_method);

    samples = sampleFunction(func, Nsamples, dvh1_distrib, dvh2_distrib,
                             alpha_distrib, beta_distrib,
                             RBEmin1_distrib, RBEmax1_distrib,
                             RBEmin2_distrib, RBEmax2_distrib);
else
    func = @(alpha, beta, RBEmin1, RBEmax1, RBEmin2, RBEmax2) CalculateMeanRR(
                dvh1list, dvh2list, opts, n1, n2, alpha, beta, RBEmin1, RBEmax1,
                RBEmin2, RBEmax2, bootstrap_samples, bootstrap_method);

    samples = sampleFunction(func, Nsamples, alpha_distrib, beta_distrib,
                             RBEmin1_distrib, RBEmax1_distrib,
                             RBEmin2_distrib, RBEmax2_distrib);
end
result = cell2mat(samples);
return;


function [low, high] = getDVHBinRanges(dvh)
% Create the low and high binning range matrix for the DVH data points.
% These can be used for the distribution sampling functions to select values
% between the low and high values.

N = length(dvh(:,1));
A = dvh(2:N,1);
B = dvh(1:N-1,1);
dx = mean(abs(A - B));
M = (A + B)*0.5;
low = zeros(size(dvh));
high = zeros(size(dvh));
low(1,1) = max(B(1) - dx, 0);
low(2:N,1) = M+eps;
low(:,2) = min(max(dvh(:,2) - 1e-5, 0), 1);
high(1:N-1,1) = M-eps;
high(N,1) = max(A(N-1) + dx, 0);
high(:,2) = min(max(dvh(:,2) + 1e-5, 0), 1);

% Adjust the limits so that the uncertainty ranges do not overlap between data points.
for k = 2:length(low(:,2))
    if low(k-1,2) <= high(k,2)
        ave = 0.5*(low(k-1,2) + high(k,2));
        ave = min(max(ave, 0), 1);
        low(k-1,2) = min(ave+eps, high(k-1,2));
        high(k,2) = max(ave-eps, low(k,2));
    end
end
return;


function result = CalculateMeanRR(dvh1list, dvh2list, opts, n1, n2, alpha, beta,
                                  RBEmin1, RBEmax1, RBEmin2, RBEmax2,
                                  bootstrap_samples, bootstrap_method)
rr = zeros(length(dvh1list), 1);
for k = 1:length(rr)
    % We have to make sure we tie the distribution end points to 1 and 0 for
    % the volume ratios, otherwise the integration becomes unstable.
    dvh1 = dvh1list{k};
    dvh1 = [0, 1; dvh1; dvh1(length(dvh1),1)+1, 0];
    dvh2 = dvh2list{k};
    dvh2 = [0, 1; dvh2; dvh2(length(dvh2),1)+1, 0];
    rr(k) = RelativeRisk('LinearQuadMultiRBE', dvh1, dvh2, opts, n1, n2,
                         alpha, beta, RBEmin1, RBEmax1, RBEmin2, RBEmax2);
end
S = bootStrap(rr, bootstrap_samples, bootstrap_method)';
M = mean(S)';
A = repmat(alpha, length(M), 1);
B = repmat(beta, length(M), 1);
Rl1 = repmat(RBEmin1, length(M), 1);
Rh1 = repmat(RBEmax1, length(M), 1);
Rl2 = repmat(RBEmin2, length(M), 1);
Rh2 = repmat(RBEmax2, length(M), 1);
result = [M, A, B, Rl1, Rh1, Rl2, Rh2];
fprintf(stdout,
        '#%d\trelative risk = %g\talpha = %g\tbeta = %g\tRBEmin1 = %g\tRBEmax1 = %g\tRBEmin2 = %g\tRBEmax2 = %g\n',
        length(M), result(1,1), result(1,2), result(1,3), result(1,4), result(1,5), result(1,6), result(1,7));
fflush(stdout);
return;
