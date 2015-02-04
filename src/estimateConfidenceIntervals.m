function result = estimateConfidenceIntervals(filenames, organs, models, options)
%result = estimateConfidenceIntervals(filenames, organs, models, options)
%
% Estimates confidence intervals for different organs and models for each given
% file using quantile estimation. In addition, a number of control parameters
% such as the minimum, maximum, median, mean and standard deviation of the
% distributions are also estimated to be able to cross check the quality of the
% resultant interval.
%
%Parameters:
%   filenames - A cell array of filenames to process. If only one file should
%               be processed then one can simply pass a string instead.
%   organs - A cell array list of organ names for which to calculate intervals.
%               If only one organ is given then this can simply be a string.
%   models - A cell array list of models for which to calculate intervals. If
%               only one model is given then this can simply be a string.
%   options - A structure containing the following fields:
%       'nsamples' - The number of samples to use for the estimation.
%               (default = 1000)
%       'confidence' - The confidence level to estimate, e.g. 0.95 for a 95%
%               confidence interval. (default = 0.95)
%       'print_progress' - Integer value to indicate if progress information
%               should be printed. Higher values generate more verbosity.
%               (default = 1).
%       'save_histos' - Boolean value indicating if the distribution histograms
%               should be saved to file as EPS images. (default = 0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    Particle Therapy Project Bergen (PTPB) - tools and models for research in
%    cancer therapy using particle beams.
%
%    Copyright (C) 2014 Particle Therapy Group Bergen
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
    help estimateConfidenceIntervals;
    return;
end

if ~ exist('filenames')
  error('No filename or names provided.');
end
if ischar(filenames)
  filenames = {filenames};
end

if exist('organs')
  if ischar(organs)
    organs = {organs};
  end
else
  organs = {};
end

if exist('models')
  if ischar(models)
    models = {models};
  end
else
  models = {};
end

% Set default values.
print_progress = 0;
save_histos = 0;
nsamples = 1000;
confidence = 0.95;

% Decode the options structure.
if exist('options')
  if isfield(options, 'print_progress')
    print_progress = options.print_progress;
  end
  if isfield(options, 'nsamples')
    nsamples = options.nsamples;
  end
  if isfield(options, 'confidence')
    confidence = options.confidence;
  end
  if isfield(options, 'save_histos')
    save_histos = options.save_histos;
  end
end

lower_quantile = (1 - confidence) / 2;
upper_quantile = 1 - lower_quantile;

result = struct;

for a = 1:length(filenames)
  filename = filenames{a};
  if print_progress > 2
    printf('Trying file "%s".\n', filename);
    fflush(stdout);
  end
  data = load(filename, 'OED_samples');

  all_organs = fieldnames(data.OED_samples);
  for b = 1:length(all_organs)
    organ = all_organs{b};
    if print_progress > 2
      printf('Trying organ "%s".\n', organ);
      fflush(stdout);
    end

    % Filter out only the selected organs.
    if length(organs) > 0
      notfound = 1;
      for k = 1:length(organs)
        if strcmp(organ, organs{k})
          notfound = 0;
          break;
        end
      end
      if notfound
        continue;
      end
    end

    all_models = fieldnames(data.OED_samples.(organ));
    for c = 1:length(all_models)
      model = all_models{c};
      if print_progress > 2
        printf('Trying model "%s".\n', model);
        fflush(stdout);
      end

      % Filter out only the selected organs.
      if length(models) > 0
        notfound = 1;
        for k = 1:length(models)
          if strcmp(model, models{k})
            notfound = 0;
            break;
          end
        end
        if notfound
          continue;
        end
      end

      if print_progress > 0
        printf('Processing file = "%s", organ = "%s", model = "%s"\n', filename, organ, model);
        fflush(stdout);
      end

      x = data.OED_samples.(organ).(model);

      if save_histos
        clf;
        hist(x, ceil(sqrt(length(x))));
        [dir, fname, ext] = fileparts(filename);
        plotname = sprintf('%s_%s_%s.eps', fname, organ, model);
        print(plotname, '-deps');
      end

      % Calculate and fill the statistics.
      [r, samples] = estimateStats(x, nsamples, lower_quantile, upper_quantile);
      result.(organ).(model) = r;

      if print_progress > 1
        printf('Interval = [%.9f .. %.9f]\n', r.lower_quantile.value, r.upper_quantile.value);
        fflush(stdout);
      end

    end
  end
end
