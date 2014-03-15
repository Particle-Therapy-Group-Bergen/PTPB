function listDoses(inputArg)
%listDoses
%listDoses(data)
%listDoses(filename)
%
% Lists all organ doses as indicated in the DVH data.
% These would be nominal doses as given by the treatment planning system.
%
% If no arguments are given then DVH data must be already loaded into the
% workspace with the variable name 'DVH_data'.
%
% Alternatively the DVH data can be passed to this function as an argument.
%
% If a character string is given instead, it is assumed to be the name of a
% Matlab workspace file containing a variable called 'DVH_data', which contains
% the DVH data. The file is loaded and all organs listed.
%
% Running example: listDoses('DVHinput') or listDoses('DVHinput.mat')
%

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

if nargin == 0
    % Try fetch data from workspace if no arguments given:
    try
        data = evalin('base', 'DVH_data');
    catch
        error('Variable "DVH_data" was not found in the workspace.');
    end
    if isfield(data, 'structures')
        organs = data.structures;
    else
        error('The global variable "DVH_data" does not appear to actually contain DVH data.');
    end
else
    if ischar(inputArg)
        try
            data = load(inputArg, 'DVH_data');
        catch
            if isempty(strfind(lasterror.message, 'empty name keyword or no data found in file'))
                rethrow(lasterror);
            end
            error('The file "%s" does not contain the DVH_data variable with DVH data.', inputArg);
        end
        if isfield(data.DVH_data, 'structures')
            organs = data.DVH_data.structures;
        else
            error('The file "%s" does not appear to contain valid DVH data.', inputArg);
        end
    elseif isfield(inputArg, 'structures')
        organs = inputArg.structures;
    else
        error('The given input does not appear to be a valid DVH data.');
    end
end

printf('No.\tMin\tMean\tMedian\tMax\tOrgan Name\n');
for n = 1:length(organs)
    printf('%d\t%g\t%g\t%g\t%g\t''%s''\n', n, organs{n}.minDose,
           organs{n}.meanDose, organs{n}.medianDose, organs{n}.maxDose,
           organs{n}.structName);
end
