function result = processPatients(dvhfiles, params, organs, models, organ_name_map, N)
%result = processPatients(dvhfiles [, params, organs, models, organ_name_map, N])
%
% This is the top level entry point for performing standard model calculations
% on one or more patient files. The output will be returned as a structure
% containing the results and sampled uncertainty distributions.
%
%Parameters:
%
% dvhfiles - A cell array of file paths containing DVH data converted to Matlab
%            format (.mat) with the convertDVHtoMatlabFile.sh script.
%
% params - A structure containing various parameters to use for the calculation.
%
% organs - A cell array of organ name strings for which to make calculations.
%
% models - A cell array of response models to use in the calculations.
%
% organ_name_map - The organ name remapping table passed to the function
%                  getDoseVolumeHistogram. Refer to that function for details
%                  about this parameters structure.
%
% N - The number of samples to produce.

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
    help processPatients;
    return;
end

if ~ exist('dvhfiles')
    error('No file names given for "dvhfiles" cell array.');
end
if ~ exist('params') || length(params) == 0
    params = setupDefaultParams();
end
if ~ exist('organs') || length(organs) == 0
    organs = {};
end
if ~ exist('models') || length(models) == 0
    models = {'LNT', 'PlateauHall', 'LinExp'};
end
if ~ exist('organ_name_map')
    organ_name_map = {};
end
if ~ exist('N')
    N = 100;
end

% Check that we have the necessary fields in the params structure.
fields = {};
fields{1} = 'dose_binning_uncertainty_model';
fields{2} = 'dose_binning_uncertainty';
fields{3} = 'volume_ratio_uncertainty_model';
fields{4} = 'volume_ratio_uncertainty';
fields{5} = 'integration_methods';
fields{6} = 'interpolation_methods';
fields{7} = 'bootstrap_max_samples';
fields{8} = 'bootstrap_sample_mode';
fields{9} = 'organs';
for n = 1:length(fields)
    if ~ isfield(params, fields{n})
        error('Missing field "%s" in the params structure.', fields{n});
    end
end

integrationMethods = params.integration_methods;
interpolationMethods = params.interpolation_methods;

% Using structs as a type of dictionary.
organs_found = struct;
models_found = struct;

% Perform calculations such as computing OEDs per input file (i.e. per patient DVH).
result = struct;
for n = 1:length(dvhfiles)
    filename = dvhfiles{n};
    result.OED_per_patient{n} = struct('filename', filename);
    dvhs = getDoseVolumeHistogram(filename, organ_name_map, params, organs{:});
    organ_names = fieldnames(dvhs);
    if length(organs) > 0
        organ_names = organs;
    end
    for m = 1:length(organ_names)
        organ = organ_names{m};
        dvh = dvhs.(organ);
        dvh.dose_binning_uncertainty_model = params.dose_binning_uncertainty_model;
        dvh.volume_ratio_uncertainty_model = params.volume_ratio_uncertainty_model;
        if ~ isfield(params.organs, organ)
            warning('Missing organ "%s" in input parameters "params.organs".', organ);
            continue;
        end
        organParams = params.organs.(organ);
        for k = 1:length(models)
            model = models{k};
            if ~ isfield(organParams, model)
                warning('Missing model "%s" in input parameters "params.organs.%s".',
                        model, organ);
                continue;
            end
            modelParams = organParams.(model);
            if ~ iscell(modelParams)
                error('Expected a cell array of model parameters for "params.organs.%s.%s".',
                      organ, model);
            end
            oeds = sampleOED(dvh, N, model, integrationMethods,
                             interpolationMethods, modelParams{:});
            result.OED_per_patient{n}.organs.(organ).(model) = oeds;
            organs_found.(organ) = 1;
            models_found.(model) = 1;
        end
    end
end

max_samples = params.bootstrap_max_samples;
sample_mode = params.bootstrap_sample_mode;

% Now apply boot-strapping to organs/models across all patients.
organ_names = fieldnames(organs_found);
model_names = fieldnames(models_found);
for n = 1:length(organ_names)
    organ = organ_names{n};
    for m = 1:length(model_names)
        model = model_names{m};
        % Collect patient data into columns.
        data = [];
        for k = 1:length(dvhfiles)
            if ~ isfield(result.OED_per_patient{k}.organs, organ)
                continue
            end
            if ~ isfield(result.OED_per_patient{k}.organs.(organ), model)
                continue
            end
            data = [data result.OED_per_patient{k}.organs.(organ).(model)];
        end
        % Boot-strap each set of samples from patients (each row in data) and
        % collect all boot-strapped samples.
        [nr, nc] = size(data);
        samples = [];
        stot = 1;
        for k = 1:nr
            s = bootStrap(data(k,:), max_samples, sample_mode)';
            [snr, snc] = size(s);
            samples(1:snr,stot:stot+snc-1) = s;
            stot += snc;
        end
        result.OED_bootstrap_samples.(organ).(model) = samples;
    end
end

return;


function params = setupDefaultParams()
% Setup a set of default parameters to be used by this script.

trapzMethod = struct('name', 'trapz', 'tolerance', 1e-6);
quadMethod = struct('name', 'quad', 'tolerance', sqrt(eps));
quadvMethod = struct('name', 'quadv', 'tolerance', 1e-6);
quadlMethod = struct('name', 'quadl', 'tolerance', eps);
quadgkMethod = struct('name', 'quadgk', 'tolerance', 1e-10);
integrationMethods = {trapzMethod, quadMethod, quadvMethod, quadlMethod, quadgkMethod};
interpolationMethods = {'linear', 'pchip', 'cubic', 'spline'};

organParams = struct;
p = getParameters('plateau_threshold', 'linexp_alphas');
names = fieldnames(p);
for n = 1:length(names)
    name = names{n};
    organParams.(name).LNT = {};  % Create default entry for LNT.
    if isfield(p.(name), 'plateau_threshold')
        A = p.(name).plateau_threshold.range_low;
        C = p.(name).plateau_threshold.value;
        B = p.(name).plateau_threshold.range_high;
        model = struct('uncertainty_model', 'triangle',
                       'params', {{A, C, B}});
        organParams.(name).PlateauHall = {model};
    end
    if isfield(p.(name), 'linexp_alpha')
        A = p.(name).linexp_alpha.range_low;
        C = p.(name).linexp_alpha.value;
        B = p.(name).linexp_alpha.range_high;
        model = struct('uncertainty_model', 'triangle',
                       'params', {{A, C, B}});
        organParams.(name).LinExp = {model};
    end
end

p = getParameters('histogram_uncertainty');
params = struct('dose_binning_uncertainty_model', 'box',
                'dose_binning_uncertainty', p.dose_binning_uncertainty,
                'volume_ratio_uncertainty_model', 'box',
                'volume_ratio_uncertainty', p.volume_ratio_uncertainty,
                'integration_methods', {integrationMethods},
                'interpolation_methods', {interpolationMethods},
                'bootstrap_max_samples', 6435,
                'bootstrap_sample_mode', 'adaptive',
                'organs', organParams);
return;
