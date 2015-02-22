function result = mergePatientResults(varargin)
%result = mergePatientResults(results1, results2 [, results3 ...])
%
% Merges the output structures produced by the processPatients.m script.
%
%Parameters:
%
% resultsN - One or more results vectors produced by a call to processPatients().

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
    help mergePatientResults;
    return;
end

for n = 1:nargin
    if ~ isstruct(varargin{n})
        error('Expected a structure for parameter %d.', n);
    end
end

% Go through all the parameters given to this function and merge them together.
result = varargin{1};
for n = 2:nargin
    results2 = varargin{n};
    fields = fieldnames(results2);
    for n = 1:length(fields)
        field = fields{n};
        if iscell(results2.(field))
            for m = 1:length(results2.(field))
                k = findSamePatient(result.(field), results2.(field){m}.filename);
                X = mergeOrganData(result.(field){k}, results2.(field){m});
                result.(field){k} = X;
            end
        else
            result.(field);
        end
    end
end
return;


function k = findSamePatient(data, filename)
% Find the first structure in the cell array of structures 'data' that has the
% same 'filename' field as the filename string given.
for n = 1:length(data)
    if strcmp(data{n}.filename, filename)
        k = n;
        return;
    end
end
% Default to a new element at the end of the array if nothing found.
k = length(data)+1;
return;


function y = mergeOrganData(a, b)
% Merge organ/model sample data from structures 'a' and 'b' into the result 'y'.
y = a;
y.filename = b.filename;
organs = fieldnames(b.organs);
for n = 1:length(organs)
    organ = organs{n};
    models = fieldnames(b.organs.(organ));
    for m = 1:length(models)
        model = models{m};
        k = length(y.organs.(organ).(model)) + 1;
        N = length(b.organs.(organ).(model));
        y.organs.(organ).(model)(k:k+N-1) = b.organs.(organ).(model);
    end
end
return;
