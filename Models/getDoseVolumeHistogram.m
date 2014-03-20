function result = getDoseVolumeHistogram(filename, varargin)
%result = getDoseVolumeDistribution(filename, [organ, ...])
%result = getDoseVolumeDistribution(filename, name_map, [organ, ...])
%
%Retrieve the dose volume histograms for a given set of organs from the
%converted DVH data file. The file should be a .mat file converted with
%the convertDVHtoMatlabFile.sh tool.
%The result will contain the dose volume histogram data points including
%their uncertainty ranges as [range_low .. range_high]. These will be stored
%as a structure with each field naming a specific organ.
%The data points are returned as a Nx2 matrix with the first column containing
%the dose values and the second the volume fraction.
%
%Parameters:
% filename - The name of the .mat file to load the data from.
% name_map - A Nx2 cell array forming key-value pairs that map the names of the
%            organ structures in the file to standard values. The first column
%            of the cell array is treated as the key (the names in the file) and
%            the second column as the values (the standard names to map to).
% organ - One or more optional organ names can be passed to indicate specific
%         structures to load from the file. If no organs are explicitly given
%         then all are loaded and returned.
%

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

% Authors: Artur Szostak <artursz@iafrica.com>

% Get the organ name map if it was given in the list of arguments.
if length(varargin) > 0 && iscell(varargin{1})
  poffset = 1;
  organ_name_map = varargin{1};
  if size(organ_name_map, 2) ~= 2
    error('The name_map parameter must be a Nx2 cell matrix.');
  end
  search_names = varargin(2:length(varargin));
else
  poffset = 0;
  search_names = varargin;
end

% Create a structure of organ names to find, and initially mark them as not found.
% Will update this later to 1 as the organs are found in the file.
for n = 1:length(search_names)
  if ~ ischar(search_names{n})
    error('Expected parameter %d to be a string; an organ name to search for in the data file.', n+poffset);
  end
  search_list.(search_names{n}) = 0;
end

% Load the data from file.
if ~ exist('filename')
  error('No file name given to load data from. Expect a converted .mat file with DVH data.');
end
data = load(filename, 'DVH_data');

params = getParameters('histogram_uncertainty');

% Extract the dose volume histograms:
result = struct;
organs = data.DVH_data.structures;
for n = 1:length(organs)
  organ_name = organs{n}.structName;
  extra_error = 0;

  % Try map the name to a standard one if we have a mapping structure.
  if exist('organ_name_map')
    pos = strmatch(organ_name, {organ_name_map{:,1}}, 'exact');
    if length(pos) > 0
      organ_name = organ_name_map{pos(1),2};
    end
  end

  % If we have a search list of organs to find, then check if this current one is in it.
  % If not then go to the next one.
  if exist('search_list')
    if isfield(search_list, organ_name)
      search_list.(organ_name) = 1;  % Mark as found.
    else
      continue;
    end
  end

  % Perform sanity checks and prepare the volume ratio from structureVolume if missing ratioToTotalVolume.
  if ~ isfield(organs{n}, 'dose')
    warning('The dose field could not be found so "%s" will be skipped.', organ_name);
    continue;
  end
  if ~ isfield(organs{n}, 'ratioToTotalVolume')
    warning('The ratioToTotalVolume field could not be found so "%s" will be skipped.', organ_name);
    continue;
  end
  if length(organs{n}.dose) == 0
    warning('The dose field is empty so "%s" will be skipped.', organ_name);
    continue;
  end
  if length(organs{n}.ratioToTotalVolume) == 0
    if ~ isfield(organs{n}, 'structureVolume')
      warning('The ratioToTotalVolume field is empty so "%s" will be skipped.', organ_name);
      continue;
    end
    warning('The ratioToTotalVolume field is empty so will use structureVolume for "%s" instead.', organ_name);
    extra_error = abs(organs{n}.structureVolume(1) - organs{n}.volume) / min([organs{n}.structureVolume(1), organs{n}.volume]);
    organs{n}.ratioToTotalVolume = organs{n}.structureVolume / organs{n}.structureVolume(1);
  end
  if length(size(organs{n}.dose)) ~= length(size(organs{n}.ratioToTotalVolume))
    warning('The dose and ratioToTotalVolume fields are not the same size so "%s" will be skipped.', organ_name);
    continue;
  end
  if size(organs{n}.dose) ~= size(organs{n}.ratioToTotalVolume)
    warning('The dose and ratioToTotalVolume fields are not the same size so "%s" will be skipped.', organs{n}.structName);
    continue;
  end

  % Rescale the volume fraction to ratio if set as percent.
  if max(organs{n}.ratioToTotalVolume) == 100
    organs{n}.ratioToTotalVolume = organs{n}.ratioToTotalVolume / 100;
  end
  % Perform sanity checks on the data point ranges:
  if min(organs{n}.ratioToTotalVolume) < 0 || 1 < max(organs{n}.ratioToTotalVolume)
    warning('Found volume ratio data points that are outside the valid range for organ "%s".', organ_name);
  end
  if min(organs{n}.dose) < 0
    warning('Found dose data points that are outside the valid range for organ "%s".', organ_name);
  end

  result.(organ_name).datapoints = [organs{n}.dose; organs{n}.ratioToTotalVolume];

  % Calculate uncertainty ranges for the data points.
  dose_min = organs{n}.dose - params.dose_binning_uncertainty;
  dose_min = dose_min .* (dose_min > 0);   % force lower bound to zero.
  dose_max = organs{n}.dose + params.dose_binning_uncertainty;
  ratio_min = organs{n}.ratioToTotalVolume - params.volume_ratio_uncertainty - extra_error;
  ratio_min = ratio_min .* (0 < ratio_min);   % force lower bound to zero.
  ratio_max = organs{n}.ratioToTotalVolume + params.volume_ratio_uncertainty + extra_error;
  ratio_max = ratio_max .* (ratio_max < 1) + (ratio_max >= 1);   % force upper bound to 1.

  % Adjust the limits so that the uncertainty ranges do not overlap between data points.
  for k = 2:length(organs{n}.dose)
    if dose_max(k-1) > dose_min(k)
      ave = 0.5*(dose_max(k-1) + dose_min(k));
      if ave < organs{n}.dose(k-1)
        ave = organs{n}.dose(k-1);
      end
      if ave > organs{n}.dose(k)
        ave = organs{n}.dose(k);
      end
      dose_max(k-1) = ave;
      dose_min(k) = ave;
    end
    if ratio_min(k-1) < ratio_max(k)
      ave = 0.5*(ratio_min(k-1) + ratio_max(k));
      if ave > organs{n}.ratioToTotalVolume(k-1)
        ave = organs{n}.ratioToTotalVolume(k-1);
      end
      if ave < organs{n}.ratioToTotalVolume(k)
        ave = organs{n}.ratioToTotalVolume(k);
      end
      ratio_min(k-1) = ave;
      ratio_max(k) = ave;
    end
  end

  result.(organ_name).range_low = [dose_min; ratio_min];
  result.(organ_name).range_high = [dose_max; ratio_max];

  % Transpose datapoints if needed to get a Nx2 matrix.
  if size(result.(organ_name).datapoints)(2) ~= 2
    result.(organ_name).datapoints = result.(organ_name).datapoints';
  end
  if size(result.(organ_name).range_low)(2) ~= 2
    result.(organ_name).range_low = result.(organ_name).range_low';
  end
  if size(result.(organ_name).range_high)(2) ~= 2
    result.(organ_name).range_high = result.(organ_name).range_high';
  end
end

% Check if any organs we were looking for were not found and print a warning about that.
if exist('search_list')
  names = fieldnames(search_list);
  for n = 1:length(names)
    if search_list.(names{n}) == 0
      warning('The organ "%s" was not found in the file "%s".', names{n}, filename);
    end
  end
end
