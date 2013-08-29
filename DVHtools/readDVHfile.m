function data = readDVHfile(filename)
%data = readDVHfile(filename)
%
%Reads in a standard DVH output text file given by the filename and returns the
%contents in the header and structures cell array.
%
%data - this contains a structure of the parsed DVH file. It will contain two
%       fields, 'header' and 'structures', structured as follows:
%  header - this contains a Matlab struct of field values that form the global
%           header of the DVH file.
%  structures - this is a cell array of Matlab structs that each contain fields
%               relevant to a biological structure (organ).

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


data = {};

% Open the file, read it in line by line and close it again.
[file, message] = fopen(filename, 'r');
if file == -1
    error('Failed to read DVH file "%s": %s', filename, message);
end
lines = {};
lineno = 0;
while ! feof(file)
    line = fgetl(file);
    lineno++;
    if ferror(file)
        error('Failed to read DVH file "%s" on line %d: %s', filename, lineno, message);
    end
    lines{lineno} = line;
end
fclose(file);

% Now parse the DVH global header
[header, linesDone] = parseTopHeader(lines);
data = struct('header', header, 'structures', {{}});

% Parse the structure data for each of the organs.
startLine = linesDone + 1;
structNum = 1;
while startLine < length(lines)
    [structData, linesDone] = parseStructure(lines, startLine);
    if length(structData) == 0
        warning('Only part of a structure was parsed. Current line %d.', linesDone);
        break;
    end
    data.structures{structNum} = structData;
    structNum++;
    startLine = linesDone + 1;
end

return;


function y = convertToNumber(x)
%Converts a string to a number taking care of cases where the string is set to 'N/A',
%in which case NaN is returned.

if strcmp(x, 'N/A')  % returns 1 if equal
    y = nan;
else
    y = str2num(x);
end
return;


function [header, linesDone] = parseTopHeader(lines)
%Parse the top level header from the input lines.
%lines should be a cellarray of all the lines read from the file.
%header will contain a structure of the header read and linesDone will indicate
%how many lines were processed by this routine.

% If we dont have any lines to process then exit now.
if length(lines) == 0
    header = {};
    linesDone = 0;
    return;
end

% Now start parsing the header keywords.
patientName = '';
patientId = '';
comment = '';
progDate = '';
exportuser = '';
progType = '';
description = '';
plan = '';
prescribedDose = nan;
prescribedDoseUnit = '';
percentForDose = nan;
percentForDoseUnit = '';
planStatus = '';

n = 1;
while n <= length(lines)
    line = lines{n};
    n++;
    if length(line) == 0
        % Skip blank lines
        continue;
    end
    tokens = regexp(line, '^\s*Patient Name\s*:\s*(.*)\s*$', 'tokens');
    if length(tokens) > 0
        patientName = tokens{1}{1};
        continue;
    end
    tokens = regexp(line, '^\s*Patient ID\s*:\s*(.*)\s*$', 'tokens');
    if length(tokens) > 0
        patientId = tokens{1}{1};
        continue;
    end
    tokens = regexp(line, '^\s*Comment\s*:\s*(.*)\s*$', 'tokens');
    if length(tokens) > 0
        comment = tokens{1}{1};
        continue;
    end
    tokens = regexp(line, '^\s*Date\s*:\s*(.*)\s*$', 'tokens');
    if length(tokens) > 0
        progDate = tokens{1}{1};
        continue;
    end
    tokens = regexp(line, '^\s*Exported by\s*:\s*(.*)\s*$', 'tokens');
    if length(tokens) > 0
        exportuser = tokens{1}{1};
        continue;
    end
    tokens = regexp(line, '^\s*Type\s*:\s*(.*)\s*$', 'tokens');
    if length(tokens) > 0
        progType = tokens{1}{1};
        continue;
    end
    tokens = regexp(line, '^\s*Description\s*:\s*(.*)\s*$', 'tokens');
    if length(tokens) > 0
        description = tokens{1}{1};
        pos = regexp(line, ':', 'start');
        pattern = sprintf('\\s{%d,}(.*)\\s*', pos);
        tokens = regexp(lines{n}, pattern, 'tokens');
        while length(tokens) > 0
            description = [description ' ' tokens{1}{1}];
            n++;
            tokens = regexp(lines{n}, pattern, 'tokens');
        end
        continue;
    end
    tokens = regexp(line, '^\s*Plan\s*:\s*(.*)\s*$', 'tokens');
    if length(tokens) > 0
        plan = tokens{1}{1};
        continue;
    end
    tokens = regexp(line, '^\s*Prescribed dose\s*\[(.*)\]\s*:\s*(.*)\s*$', 'tokens');
    if length(tokens) > 0
        prescribedDoseUnit = tokens{1}{1};
        prescribedDose = convertToNumber(tokens{1}{2});
        continue;
    end
    tokens = regexp(line, '^\s*% for dose\s*\((.*)\)\s*:\s*(.*)\s*$', 'tokens');
    if length(tokens) > 0
        percentForDoseUnit = tokens{1}{1};
        percentForDose = convertToNumber(tokens{1}{2});
        continue;
    end
    tokens = regexp(line, '^\s*Plan Status\s*:\s*(.*)\s*$', 'tokens');
    if length(tokens) > 0
        planStatus = tokens{1}{1};
        continue;
    end
    % If we did not match and continue from any of the keywords above and
    % ended up at this point then we assume the section we are looking at
    % Is no longer the header and we stop parsing.
    n--;  % Decrement n since this last line was not processed.
    break;
end

% Check that we parsed all relevant information for the header.
% If not then warn the user.
if length(patientName) == 0 || length(patientId) == 0 || length(progDate) == 0 \
        || length(progType) == 0 || length(plan) == 0 \
        || prescribedDose == nan || length(prescribedDoseUnit) == 0 \
        || percentForDose == nan || length(percentForDoseUnit) == 0
    warning('Could not find all the relevant DVH header fields. Current line %d.', n);
end

header = struct('patientName', patientName, 'patientId', patientId,
                'comment', comment, 'progDate', progDate, 'exportuser', exportuser,
		'progType', progType, 'planStatus', planStatus,
                'description', description, 'plan', plan,
                'prescribedDose', prescribedDose,
                'prescribedDoseUnit', prescribedDoseUnit,
                'percentForDose', percentForDose,
                'percentForDoseUnit', percentForDoseUnit);
linesDone = n-1;
return;


function [data, linesDone] = parseStructure(lines, startLine)
%Parse a structure header and data from the input lines.
%lines should be a cellarray of all the lines read from the file.
%startLine should indicate the first line from which to start processing.
%data will contain the structure of the parsed header and data lines.
%linesDone will indicate how many lines were processed by this routine,
%including the skipped lines, i.e. the line on which successful parsing ended.

% If we dont have any lines to process then exit now.
if startLine > length(lines)
    data = []
    linesDone = 0;
    return;
end

% Now start parsing the header keywords.
structName = '';
approvalStatus = '';
plan = '';
course = '';
volume = nan;
volumeUnit = '';
doseCoverage = nan;
doseCoverageUnit = '';
samplingCoverage = nan;
samplingCoverageUnit = '';
minDose = nan;
minDoseUnit = '';
maxDose = nan;
maxDoseUnit = '';
meanDose = nan;
meanDoseUnit = '';
modalDose = nan;
modalDoseUnit = '';
medianDose = nan;
medianDoseUnit = '';
standardDeviation = nan;
standardDeviationUnit = '';
equivSphereDiam = nan;
equivSphereDiamUnit = '';
conformityIndex = '';
gradientMeasure = nan;
gradientMeasureUnit = '';

n = startLine;
while n <= length(lines)
    line = lines{n};
    n++;
    if length(line) == 0
        % Skip blank lines
        continue;
    end
    tokens = regexp(line, '^\s*Structure\s*:\s*(.*)\s*$', 'tokens');
    if length(tokens) > 0
        structName = tokens{1}{1};
        continue;
    end
    tokens = regexp(line, '^\s*Approval Status\s*:\s*(.*)\s*$', 'tokens');
    if length(tokens) > 0
        approvalStatus = tokens{1}{1};
        continue;
    end
    tokens = regexp(line, '^\s*Plan\s*:\s*(.*)\s*$', 'tokens');
    if length(tokens) > 0
        plan = tokens{1}{1};
        continue;
    end
    tokens = regexp(line, '^\s*Course\s*:\s*(.*)\s*$', 'tokens');
    if length(tokens) > 0
        course = tokens{1}{1};
        continue;
    end
    tokens = regexp(line, '^\s*Volume\s*\[(.*)\]\s*:\s*(.*)\s*$', 'tokens');
    if length(tokens) > 0
        volumeUnit = tokens{1}{1};
        volume = convertToNumber(tokens{1}{2});
        continue;
    end
    tokens = regexp(line, '^\s*Dose Cover.\s*\[(.*)\]\s*:\s*(.*)\s*$', 'tokens');
    if length(tokens) > 0
        doseCoverageUnit = tokens{1}{1};
        doseCoverage = convertToNumber(tokens{1}{2});
        continue;
    end
    tokens = regexp(line, '^\s*Sampling Cover.\s*\[(.*)\]\s*:\s*(.*)\s*$', 'tokens');
    if length(tokens) > 0
        samplingCoverageUnit = tokens{1}{1};
        samplingCoverage = convertToNumber(tokens{1}{2});
        continue;
    end
    tokens = regexp(line, '^\s*Min Dose\s*\[(.*)\]\s*:\s*(.*)\s*$', 'tokens');
    if length(tokens) > 0
        minDoseUnit = tokens{1}{1};
        minDose = convertToNumber(tokens{1}{2});
        continue;
    end
    tokens = regexp(line, '^\s*Max Dose\s*\[(.*)\]\s*:\s*(.*)\s*$', 'tokens');
    if length(tokens) > 0
        maxDoseUnit = tokens{1}{1};
        maxDose = convertToNumber(tokens{1}{2});
        continue;
    end
    tokens = regexp(line, '^\s*Mean Dose\s*\[(.*)\]\s*:\s*(.*)\s*$', 'tokens');
    if length(tokens) > 0
        meanDoseUnit = tokens{1}{1};
        meanDose = convertToNumber(tokens{1}{2});
        continue;
    end
    tokens = regexp(line, '^\s*Modal Dose\s*\[(.*)\]\s*:\s*(.*)\s*$', 'tokens');
    if length(tokens) > 0
        modalDoseUnit = tokens{1}{1};
        modalDose = convertToNumber(tokens{1}{2});
        continue;
    end
    tokens = regexp(line, '^\s*Median Dose\s*\[(.*)\]\s*:\s*(.*)\s*$', 'tokens');
    if length(tokens) > 0
        medianDoseUnit = tokens{1}{1};
        medianDose = convertToNumber(tokens{1}{2});
        continue;
    end
    tokens = regexp(line, '^\s*STD\s*\[(.*)\]\s*:\s*(.*)\s*$', 'tokens');
    if length(tokens) > 0
        standardDeviationUnit = tokens{1}{1};
        standardDeviation = convertToNumber(tokens{1}{2});
        continue;
    end
    tokens = regexp(line, '^\s*Equiv. Sphere Diam.\s*\[(.*)\]\s*:\s*(.*)\s*$', 'tokens');
    if length(tokens) > 0
        equivSphereDiamUnit = tokens{1}{1};
        equivSphereDiam = convertToNumber(tokens{1}{2});
        continue;
    end
    tokens = regexp(line, '^\s*Conformity Index\s*:\s*(.*)\s*$', 'tokens');
    if length(tokens) > 0
        conformityIndex = tokens{1}{1};
        continue;
    end
    tokens = regexp(line, '^\s*Gradient Measure\s*\[(.*)\]\s*:\s*(.*)\s*$', 'tokens');
    if length(tokens) > 0
        gradientMeasureUnit = tokens{1}{1};
        gradientMeasure = convertToNumber(tokens{1}{2});
        continue;
    end
    % If we did not match and continue from any of the keywords above and
    % ended up at this point then we assume the section we are looking at
    % Is no longer the header and we stop parsing.
    n--;  % Decrement n since this last line was not processed.
    break;
end

if length(structName) == 0
    data = {};
    linesDone = n-1;
    return;
end

% What is left to process is the data lines which is a table of numbers.
dose = [];
doseUnit = '';
relativeDose = [];
relativeDoseUnit = '';
ratioToTotalVolume = [];
ratioToTotalVolumeUnit = '';
row = 1;

while n <= length(lines)
    line = lines{n};
    n++;
    if length(line) == 0
        % Skip blank lines
        continue;
    end
    tokens = regexp(line, '^\s*Dose\s*\[(.*)\]\s*Relative dose\s*\[(.*)\]\s*Ratio of Total Structure Volume\s*\[(.*)\]\s*$', 'tokens');
    if length(tokens) > 0
        doseUnit = tokens{1}{1};
        relativeDoseUnit = tokens{1}{2};
        ratioToTotalVolumeUnit = tokens{1}{3};
        continue;
    end
    tokens = regexp(line, '^\s*([+\-\.eE0-9]+)\s+([+\-\.eE0-9]+)\s+([+\-\.eE0-9]+)\s*$', 'tokens');
    if length(tokens) > 0
        dose(row) = convertToNumber(tokens{1}{1});
        relativeDose(row) = convertToNumber(tokens{1}{2});
        ratioToTotalVolume(row) = convertToNumber(tokens{1}{3});
        row++;
        continue;
    end
    % If we did not match and continue from any of the regexp checks above
    % and ended up at this point then we assume the section we are looking
    % at is no longer part of the data table and we stop parsing.
    n--;  % Decrement n since this last line was not processed.
    break;
end

% Check that we parsed all relevant information for the structure.
% If not then warn the user.
if length(structName) == 0 || length(approvalStatus) == 0 || length(plan) == 0 \
        || length(course) == 0 || volume == nan || length(volumeUnit) == 0 \
        || doseCoverage == nan || length(doseCoverageUnit) == 0 \
        || samplingCoverage == nan || length(samplingCoverageUnit) == 0 \
        || minDose == nan || length(minDoseUnit) == 0 \
        || maxDose == nan || length(maxDoseUnit) == 0 \
        || meanDose == nan || length(meanDoseUnit) == 0 \
        || modalDose == nan || length(modalDoseUnit) == 0 \
        || medianDose == nan || length(medianDoseUnit) == 0 \
        || standardDeviation == nan || length(standardDeviationUnit) == 0 \
        || length(doseUnit) == 0 || length(relativeDoseUnit) == 0 \
        || length(ratioToTotalVolumeUnit) == 0 \
    warning('Could not find all the relevant fields for a DVH structure. Current line %d.', n);
end

data = struct('structName', structName, 'approvalStatus', approvalStatus,
              'plan', plan, 'course', course,
              'volume', volume, 'volumeUnit', volumeUnit,
              'doseCoverage', doseCoverage, 'doseCoverageUnit', doseCoverageUnit,
              'samplingCoverage', samplingCoverage, 'samplingCoverageUnit', samplingCoverageUnit,
              'minDose', minDose, 'minDoseUnit', minDoseUnit,
              'maxDose', maxDose, 'maxDoseUnit', maxDoseUnit,
              'meanDose', meanDose, 'meanDoseUnit', meanDoseUnit,
              'modalDose', modalDose, 'modalDoseUnit', modalDoseUnit,
              'medianDose', medianDose, 'medianDoseUnit', medianDoseUnit,
              'standardDeviation', standardDeviation, 'standardDeviationUnit', standardDeviationUnit,
              'equivSphereDiam', equivSphereDiam, 'equivSphereDiamUnit', equivSphereDiamUnit,
              'conformityIndex', conformityIndex,
              'gradientMeasure', gradientMeasure, 'gradientMeasureUnit', gradientMeasureUnit,
              'dose', dose, 'doseUnit', doseUnit,
              'relativeDose', relativeDose, 'relativeDoseUnit', relativeDoseUnit,
              'ratioToTotalVolume', ratioToTotalVolume, 'ratioToTotalVolumeUnit', ratioToTotalVolumeUnit);
linesDone = n-1;
return;
