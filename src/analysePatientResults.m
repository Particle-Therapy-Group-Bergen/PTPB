function table = analysePatientResults(results, make_hists, confidence)
%table = analysePatientResults(results [, make_hists, confidence])
%
% Analyses the results produced by processPatients(). Will produce histograms
% and result tables containing means, standard deviations and confidence
% intervals.
%
%Parameters:
%
% results - The data structure produced by the function processPatients().
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

if nargin == 0
    % Print help message if no arguments are given.
    help analysePatientResults;
    return;
end

if ~ exist('make_hists')
    make_hists = 1;
end
if ~ exist('confidence')
    confidence = 0.95;
end

% Calculate the lower and upper quantiles to calculate for the confidence
% interval requested.
ci_diff = (1 - confidence) * 0.5;
lower_quantile = ci_diff;
upper_quantile = 1 - ci_diff;

table = {};
row = 1;
fields = fieldnames(results);
for n = 1:length(fields)
    field = fields{n};
    if iscell(results.(field))
        for m = 1:length(results.(field))
            [organs, models, stats] = analyseOrgansAndModels(results.(field){m}.organs, lower_quantile, upper_quantile);
            for k = 1:length(organs)
                table{row,1} = field;
                table{row,2} = results.(field){m}.filename;
                table{row,3} = organs{k};
                table{row,4} = models{k};
                table{row,5} = stats{k};
                row += 1;
            end
        end
    else
        [organs, models, stats] = calcOrganAndModelAverages(results.(field), lower_quantile, upper_quantile);
        for k = 1:length(organs)
            table{row,1} = field;
            table{row,2} = '';
            table{row,3} = organs{k};
            table{row,4} = models{k};
            table{row,5} = stats{k};
            row += 1;
        end
    end
end
printTable(table);
return;


function [organs, models, stats] = analyseOrgansAndModels(data, lower_quantile, upper_quantile)
organs = {};
models = {};
stats = {};
k = 1;
organ_names = fieldnames(data);
for n = 1:length(organ_names)
    organ = organ_names{n};
    model_names = fieldnames(data.(organ));
    for m = 1:length(model_names)
        model = model_names{m};
        organs{k} = organ;
        models{k} = model;
        stats{k} = estimateStats(data.(organ).(model), lower_quantile, upper_quantile);
        k += 1;
    end
end
return;


function [organs, models, stats] = calcOrganAndModelAverages(data, lower_quantile, upper_quantile)
organs = {};
models = {};
stats = {};
k = 1;
organ_names = fieldnames(data);
for n = 1:length(organ_names)
    organ = organ_names{n};
    model_names = fieldnames(data.(organ));
    for m = 1:length(model_names)
        model = model_names{m};
        organs{k} = organ;
        models{k} = model;
        samples = mean(data.(organ).(model))';
        stats{k} = estimateStats(samples, lower_quantile, upper_quantile);
        k += 1;
    end
end
return;


function y = estimateStats(x, lower_quantile, upper_quantile)
q = quantile(x, [lower_quantile, upper_quantile], 1, 8);
y = [min(x), max(x), q', mean(x), median(x), std(x)];
return;


function printTable(table)
[nr, nc] = size(table);
maxtype = max(cellfun('length', {table{:,1}}));
maxfilename = max(cellfun('length', {table{:,2}}));
maxorgan = max(cellfun('length', {table{:,3}}));
maxmodel = max(cellfun('length', {table{:,4}}));
printf('%*s\t%*s\t%*s\t%*s\t%12s\t%12s\t%12s\t%12s\t%12s\t%12s\t%12s\n',
       maxtype, 'Type', maxfilename, 'Filename', maxorgan, 'Organ',
       maxmodel, 'Model', 'Min', 'Max', 'CI low', 'CI high', 'Mean',
       'Median', 'Std.Dev.');
for n = 1:nr
    stats = table{n,5};
    printf('%*s\t%*s\t%*s\t%*s\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n',
           maxtype, table{n,1}, maxfilename, table{n,2}, maxorgan, table{n,3},
           maxmodel, table{n,4}, stats(1), stats(2), stats(3), stats(4),
           stats(5), stats(6), stats(7));
end
return;
