function result = sampleDVH(dvhs, N, V, params)
%result = sampleDVH(dvhs [, N, V, params])
%
% This function randomly samples DVH structures and applies a boot-strap
% procedure over all DVHs given. The returned result are the mean DVH samples
% per volume ratio bin.
%
%Parameters:
%
% dvhs - A cell array of DVH structures as returned by getDoseVolumeHistogram().
%
% N - The number of times to sample each input DVH.
%
% V - The volume ratio bins (interpolation points) to sample.
%
% params - A structure containing parameters to use for the calculation:
%            dose_binning_uncertainty_model - uncertainty model for dose bin,
%                                             default 'box'.
%            volume_ratio_uncertainty_model - uncertainty model for volume bin,
%                                             default 'box'.
%            interpolation_method - interpolation method, default 'pchip'.
%            bootstrap_max_samples - maximum samples for bootstrap method,
%                                    default 6435
%            bootstrap_sample_mode - sampling mode for bootstrap method,
%                                    default 'adaptive'

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

% Author: Artur Szostak <artursz@iafrica.com>

if nargin == 0
    % Print help message if no arguments are given.
    help sampleDVH;
    return;
end

if ~ exist('params')
    params = struct;
    p = getParameters('histogram_uncertainty');
    params.dose_binning_uncertainty_model = 'box';
    params.volume_ratio_uncertainty_model = 'box';
    params.interpolation_method = 'pchip';
    params.bootstrap_max_samples = 6435;
    params.bootstrap_sample_mode = 'adaptive';
end

dose_model = params.dose_binning_uncertainty_model;
volume_model = params.volume_ratio_uncertainty_model;

% Sample the each DVH N times.
sampled_dvhs = {};
for n = 1:length(dvhs)
    dvh = dvhs{n};
    dose = dvh.datapoints(:,1);
    doseLow = dvh.range_low(:,1);
    doseHigh = dvh.range_high(:,1);
    volume = dvh.datapoints(:,2);
    volumeLow = dvh.range_low(:,2);
    volumeHigh = dvh.range_high(:,2);
    dose_samples = sampleDistribution(dose_model, N, doseLow', doseHigh');
    volume_samples = sampleDistribution(volume_model, N, volumeLow', volumeHigh');
    sampled_dvhs{n} = struct('doses', dose_samples, 'volumes', volume_samples);
end

max_samples = params.bootstrap_max_samples;
sample_mode = params.bootstrap_sample_mode;

% Now interpolate for each sampled DVH on the volume-ratio bins indicated by V.
% This forms a group of interpolated points, one per V and DVH. Then for each
% point in V, boot-strap all values in the group and calculate the means.
result = [];
for n = 1:N
    Y = zeros(length(dvhs), length(V));
    for m = 1:length(dvhs)
        dvh_samples = [sampled_dvhs{m}.doses(n,:); sampled_dvhs{m}.volumes(n,:)]';
        Y(m,:) = doseInterpolate(V, dvh_samples, params.interpolation_method);
    end
    M = [];
    for m = 1:length(V)
        samples = bootStrap(Y(:,m), max_samples, sample_mode)';
        M = [M mean(samples)'];
    end
    result = [result; M];
end
