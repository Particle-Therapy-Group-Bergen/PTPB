function result = convertEQD(dvh, n, alpha_beta, T)
%result = convertEQD(dvh, n, alpha_beta [, T])
%
% Converts the DVH histogram's dose bins to the isoeffective dose.
% The following equations is applied to each bin i:
%
%               D_i/n + alpha/beta
%  EQD_i = D_i --------------------
%                 T + alpha/beta
%
% Parameters:
%  dvh - The DVH structure as returned by getDoseVolumeHistogram(). Note that it
%        will normally return a structure for all the organs. Thus one needs to
%        select which organ to actually use. For example, if we have a input
%        file called 'mydata.mat' which was generated with the script
%        convertDVHtoMatlabFile.sh, we would use the following code to apply the
%        conversion to the Lungs:
%          data = getDoseVolumeHistogram('mydata.mat');
%          dvh = convertEQD(data.Lungs, n, alpha_beta, T);
%  n - This is the number of fractionations applied in the treatment plan.
%  alpha_beta - Is the alpha/beta ratio value.
%  T - Is the total dose for one fractionation. This parameter is optional and
%      a value of 2 Gy is used by default.
%
% Note: The input structure for the 'dvh' parameter is a struct() object with
% the following 3 fields: 'datapoints', 'range_high', 'range_low'. Each field is
% a 2 column matrix and all have identical sizes. The first column indicates the
% dose value and the second the volume fraction. Each row of 'range_high'
% indicates the upper error bar or maximal valid value for the corresponding row
% in 'datapoints'. Similarly 'range_low' indicates the lower valid value.
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

if ~ exist('T')
  T = 2;
end

datapoints = [EQD(dvh.datapoints(:,1), n, alpha_beta, T), dvh.datapoints(:,2)];
range_low = [EQD(dvh.range_low(:,1), n, alpha_beta, T), dvh.range_low(:,2)];
range_high = [EQD(dvh.range_high(:,1), n, alpha_beta, T), dvh.range_high(:,2)];
result = struct('datapoints', datapoints, 'range_low', range_low, 'range_high', range_high);
return;


function result = EQD(D, n, alpha_beta, T)
% Apply EQD calculation.
result = D.* ((D./n + alpha_beta) ./ (T + alpha_beta));
return;
