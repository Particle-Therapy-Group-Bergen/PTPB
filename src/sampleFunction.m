function [results, samples] = sampleFunction(func, varargin)
%result = sampleFunction(func, ...)
%result = sampleFunction(func, N, ...)
%
% This function will sample input distributions to get random input values to
% apply to the given function. This is performed N times and the calculated
% results are returned as a cell array of size N, one element for each set of
% sampled input values.
%
%Parameters:
% func - is the function to calculate with.
%
% N - optional value indicating the number of times to sample the input
% distributions and calculate a result. If not given then only one simple is
% calculated. Note: explicitly giving a value for N will always cause the output
% to be placed in a cell array, even for N=1.
%
% All subsequent input parameters must be structures that define the input
% distributions to sample, one structure per input parameter to 'func'. The
% structure must contain the following fields:
%   'type' - indicates the type of distribution to use. This can be one of
%            'struct', 'array', 'matrix', 'vector' or any distribution name
%            accepted by the sampleDistribution() function.
%   'params' - this is either a structure, cell array of parameters passed to
%              the sampleDistribution() function, or values to select from,
%              depending on the value 'type'.
%   'size' - a vector of integers indicating the size of the resultant matrix or
%            cell array of samples to produce. For example, if set to [2, 3] then
%            a matrix with 2 rows and 3 columns of samples is produced. Cell
%            arrays are produced instead if 'type' is 'struct' or 'array'. This
%            field is optional and assumed to be 1 (i.e. a scalar) if not given.
% The different values for 'type' are handled in the following manner:
%   'struct' - allows to produce structure objects with each field filled from
%              randomly sampled values. The 'params' field will have to be a
%              structure, where each field's value is again a structure as
%              defined above. This allows to produce structures of any depth and
%              complexity. The field names for the produced structure will be the
%              same as those of the 'params' structure.
%   'array' - will produce a cell array of randomly sampled values. In this case
%             the 'params' field must be a cell array of structures. Each one of
%             the same format as defined above, to indicate the distribution to
%             sample for that cell of the array.
%   'matrix' - This is the same concept as for 'array', except the results from
%              each distribution are concatenated into a matrix instead.
%   'matrixT' - Similar to 'matrix', but the resultant matrix is transposed.
%   'vector' - in this case the 'params' field must be a vector of values from
%              which an entry will be randomly selected. The selection uses a
%              uniform random number, thus to have higher probabilities of
%              certain values one should put in duplicates.
%   <distrib> - if a distribution name supported by sampleDistribution() is given
%               then the 'params' field must be a cell array of parameters to
%               pass to sampleDistribution().
%
% The output will have the N results of func(samples{i}{:}) placed into the
% 'result' cell array with the corresponding input values in the samples cell
% array.
%
%Example:
% TODO
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

% Authors: Artur Szostak <artursz@iafrica.com>

if nargin == 0
    % Print help message if no arguments are given.
    help sampleFunction;
    return;
end

if length(varargin) > 0 && isscalar(varargin{1}) && isnumeric(varargin{1})
    N = varargin{1};
    unpack = false;
    definitions = varargin(2:length(varargin));
    Dstart = 2;
else
    N = 1;
    unpack = true;
    definitions = varargin;
    Dstart = 1;
end

results = cell(N, 1);
samples = cell(N, 1);
for n = 1:N
    M = length(definitions);
    samplerow = cell(M, 1);
    for m = 1:M
        try
            samplerow{m} = handleDistribDefinition(definitions{m});
        catch
            err = lasterror();
            err.message = sprintf('Problem with input parameter %d: %s',
                                  Dstart+m, err.message);
            rethrow(err);
        end
    end
    samples{n} = samplerow;
    results{n} = func(samplerow{:});
end
if unpack
    results = results{1};
    samples = samples{1};
end
return;


function result = handleDistribDefinition(distrib)
if ~ isstruct(distrib)
    error('Distribution definition must be a structure created with struct().');
end
if ~ isfield(distrib, 'type')
    error('Distribution definition is missing the "type" field.');
end
if ~ isfield(distrib, 'params')
    error('Distribution definition is missing the "params" field.');
end
params = distrib.params;
sizeval = [1 1];
if isfield(distrib, 'size')
    sizeval = distrib.size;
    if ~ (isvector(sizeval) && ismatrix(sizeval))
        error('The "size" field must be a vector.');
    end
    if length(sizeval) == 1
        sizeval = [sizeval 1];
    end
end
N = prod(sizeval);
switch distrib.type
    case 'struct'
        if ~ isstruct(params)
            error('Expect a structure for "params" when "type" is "struct".');
        end
        if N > 1
            values = cell(N, 1);
            for n = 1:N
                values{n} = structfun(@handleDistribDefinition, params,
                                      'UniformOutput', false);
            end
            result = reshape(values, sizeval);
        else
            result = structfun(@handleDistribDefinition, params,
                               'UniformOutput', false);
        end
    case 'array'
        if ~ iscell(params)
            error('Expect a cell array for "params" when "type" is "array".');
        end
        if N > 1
            values = cell(N, 1);
            for n = 1:N
                values{n} = cellfun(@handleDistribDefinition, params,
                                    'UniformOutput', false);
            end
            result = reshape(values, sizeval);
        else
            result = cellfun(@handleDistribDefinition, params,
                             'UniformOutput', false);
        end
    case 'matrix'
        if ~ iscell(params)
            error('Expect a cell array for "params" when "type" is "matrix".');
        end
        if N > 1
            values = cell(N, 1);
            for n = 1:N
                values{n} = cell2mat(cellfun(@handleDistribDefinition, params,
                                             'UniformOutput', false));
            end
            result = reshape(values, sizeval);
        else
            result = cell2mat(cellfun(@handleDistribDefinition, params,
                                      'UniformOutput', false));
        end
    case 'matrixT'
        if ~ iscell(params)
            error('Expect a cell array for "params" when "type" is "matrix".');
        end
        if N > 1
            values = cell(N, 1);
            for n = 1:N
                values{n} = cell2mat(cellfun(@handleDistribDefinition, params,
                                             'UniformOutput', false))';
            end
            result = reshape(values, sizeval);
        else
            result = cell2mat(cellfun(@handleDistribDefinition, params,
                                      'UniformOutput', false))';
        end
    case 'vector'
        imax = length(params);
        index = floor(rand(N, 1) * imax) + 1;
        if iscell(params)
            values = params{index};
        else
            values = params(index);
        end
        if N > 1
            result = reshape(values, sizeval);
        else
            result = values;
        end
    otherwise
        if ~ iscell(params)
            error('Expect a cell array for "params" when "type" is "%s".',
                  distrib.type);
        end
        if N > 1
            values = sampleDistribution(distrib.type, N, params{:});
            result = reshape(values, sizeval);
        else
            result = sampleDistribution(distrib.type, N, params{:});
        end
end
return;
