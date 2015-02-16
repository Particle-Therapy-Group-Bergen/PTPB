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
if ~ exist('organ')
    error('No organ name given.');
end
if ~ exist('intervals')
    intervals = [0, 0.5, 1, 10,  13, 100];
end
if ~ exist('interpolation_method')
    interpolation_method = 'pchip';
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

dvh = getDoseVolumeHistogram(filename, organ_name_map, organ).(organ);
x = intervals;
y = InterpolateDVH(dvh, x, interpolation_method);
N = length(y);
volumeFractions = y(1:N-1) - y(2:N);
bar(volumeFractions * 100);
title(sprintf('Volume fractions in given dose range for %s', organ));
xlabel('Dose ranges [Gy]');
ranges = {};
for n = 1:N-1
    ranges{n} = sprintf('%g - %g', x(n), x(n+1));
end
set(gca, 'XTick', 1:N-1, 'XTickLabel', ranges);
ylabel('Volume fraction [%]');
result = volumeFractions;
return;


function y = InterpolateDVH(dvh, x, interpolation_method)
% Interpolation function to interpolate along the dose axis.
datapoints = dvh.datapoints';
y = interp1(datapoints(1,:), datapoints(2,:), x, interpolation_method);
y(find(isnan(y))) = 0;
return;
