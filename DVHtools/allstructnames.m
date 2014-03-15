function allstructnames(filename)
% Writes all structure names in a DVH file to screen.
%
% For example, from within MatLab/Octave console run:
%   allstructnames('InputDVHfile')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    Particle Therapy Project Bergen (PTPB) - tools and models for research in
%    cancer therapy using particle beams.
%
%    Copyright (C) 2013 Particle Therapy Group Bergen
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

% printing all organ names

% reads DVHinput and prints patient ID
fprintf('Reading file "%s"...\n', filename);
fflush(stdout);
data = readDVHfile(filename);
fprintf('Done.\n');
fflush(stdout);
fprintf('DVH for patient: %s\n', data.header.patientName);

for j=1:length(data.structures)
	structure=data.structures{j}.structName;
	fprintf('number %d \n', j);
	fprintf('structure %s \n', structure);
end


