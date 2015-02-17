function result = calcVolumeFractions(filename, organ, intervals, interpolation_method, organ_name_map)
%result = calcVolumeFractions(filename, organ [, intervals, interpolation_method, organ_name_map])
%
% This is an example script to show how to calculate the volume fractions receiving
% a doses in a certain range.
%
% filename - This must be a string with the name of a file that contains the DVH data.
%            The file must have been converted with convertDVHtoMatlabFile.sh.
%
% organ - A string with the name of the organ structure to process.
%
% intervals - A vector of dose values defining the boundary of the dose ranges.
%             For example [0, 1, 3] defines two intervals, 0 to 1 Gy and 1 to 3 Gy.
%
% interpolation_method - The interpolation method to use in the interp1 function.
%                        By default this is 'pchip'.
%
% organ_name_map - The organ name remapping table as passed to the function
%                  getDoseVolumeHistogram(). See that function for more details.
%
% Example usage:
%
%   calcVolumeFractions('3dcrt12.mat', 'Bladder')
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

% Author: Artur Szostak <artursz@iafrica.com>

if ~ exist('filename')
    error('No DVH filename given.');
end
if ~ ischar(filename)
    error('The filename must be a string.')
end
if ~ exist('organ')
    error('No organ name given.');
end
if ~ ischar(organ)
    error('The organ name must be a string.')
end
if ~ exist('intervals')
    intervals = [0, 0.5, 1, 10,  13, 100];
end
if ~ isvector(intervals) || length(intervals) < 2
    error('The intervals argument must be a vector of at least 2 elements to define a range.');
end
if ~ exist('interpolation_method')
    interpolation_method = 'pchip';
end
if ~ ischar(interpolation_method)
    error('The interpolation_method parameter must be a string.')
end
if ~ exist('organ_name_map')
    organ_name_map = {
            'BODY',               'Body';
            'Cribriform plate',   'CribriformPlate';
            'CribiformPlate',     'CribriformPlate';
            'Rectum_P_MT',        'Colon';
            'Bladder_P',          'Bladder';
        };
end
if ~ iscell(organ_name_map)
    error('The organ_name_map parameter must be a cell array.')
end

dvh = getDoseVolumeHistogram(filename, organ_name_map, organ).(organ);
x = intervals;
y = InterpolateDVH(dvh, x, interpolation_method);
N = length(y);
volumeFractions = y(1:N-1) - y(2:N);
% Draw a bar graph only if there is more than one range.
% Otherwise the bar() function fails when it gets a single scalar.
if N > 2
    bar(volumeFractions * 100);
    title(sprintf('Volume fractions in given dose range for %s', organ));
    xlabel('Dose ranges [Gy]');
    ranges = {};
    for n = 1:N-1
        ranges{n} = sprintf('%g - %g', x(n), x(n+1));
    end
    set(gca, 'XTick', 1:N-1, 'XTickLabel', ranges);
    ylabel('Volume fraction [%]');
end
result = volumeFractions;
return;


function y = InterpolateDVH(dvh, x, interpolation_method)
% Interpolation function to interpolate along the dose axis.
datapoints = dvh.datapoints';
y = interp1(datapoints(1,:), datapoints(2,:), x, interpolation_method);
y(find(isnan(y))) = 0;
return;
