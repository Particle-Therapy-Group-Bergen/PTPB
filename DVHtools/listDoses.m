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
