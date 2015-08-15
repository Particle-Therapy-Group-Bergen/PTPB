function result = applyNeutronCorrection(dvh, organ)
%result = applyNeutronCorrection(dvh, organ)
%
%Parameters:
%  dvh - The DVH structure as returned by getDoseVolumeHistogram(). Note that it
%        will normally return a structure for all the organs. Thus one needs to
%        select which organ to actually use. For example, if we have a input
%        file called 'mydata.mat' which was generated with the script
%        convertDVHtoMatlabFile.sh, we would use the following code to apply the
%        conversion to the Lungs:
%          data = getDoseVolumeHistogram('mydata.mat');
%          dvh = applyNeutronCorrection(data.Lungs, 'Lungs');
%  organ - A string indicating the name of the organ the DVH is for.
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

datapoints = [scale_dose(dvh.datapoints(:,1), organ), dvh.datapoints(:,2)];
range_low = [scale_dose(dvh.range_low(:,1), organ), dvh.range_low(:,2)];
range_high = [scale_dose(dvh.range_high(:,1), organ), dvh.range_high(:,2)];
result = struct('datapoints', datapoints, 'range_low', range_low, 'range_high', range_high);
if isfield(dvh, 'units')
    result.units = dvh.units;
end
return


function result = scale_dose(dose, organ)
% Apply Neutron dose scaling

scale_ratio = 23.4/30.6/1000*0.17;
conf_factor = 1.959963984540;
switch organ
  case 'Lungs'
    result = dose + 26.91 * scale_ratio * (1+0.25e-2 * conf_factor);
  case 'Stomach'
    result = dose + 18.41 * scale_ratio * (1+0.22e-2 * conf_factor);
  case 'Colon'
    result = dose + 13.95 * scale_ratio * (1+0.26e-2 * conf_factor);
  case 'Bladder'
    result = dose + 6.35 * scale_ratio * (1+0.67e-2 * conf_factor);
  case 'Liver'
    result = dose + 18.91 * scale_ratio * (1+0.12e-2 * conf_factor);
  case 'Thyroid'
    result = dose + 36.88 * scale_ratio * (1+0.54e-2 * conf_factor);
  case 'Breast'
    result = dose + 19.37 * scale_ratio * (1+1.13e-2 * conf_factor);
  case 'Prostate'
    result = dose + 5.79 * scale_ratio * (1+1.29e-2 * conf_factor);
  otherwise
    error('Neutron dose not known for organ %s', organ);
end

return
