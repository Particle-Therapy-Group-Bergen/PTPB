function result = getParameters(varargin)
%result = getParameters(name, [name, ...])
%
% Returns model parameters and their uncertainty ranges used in dose and risk calculations.
% Parameters are returned by named fields in a structure, including their uncertainty range
% as [range_low .. range_high].
% The one or more names of parameter sets to return can can be given, in which case only a subset
% of the parameters is returned. If no names are given then all parameters are returned.
% Possible named parameter sets include:
%   histogram_uncertainty - Parameters indicating the uncertainties for the histogram binning.
%   plateau_threshold - The threshold values for Hall's plateau model.
%   linexp_alphas - The alpha parameters for the linear-exponential model.
%   competition_params - Parameters for the competition model.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    Particle Therapy Project Bergen (PTPB) - tools and models for research in
%    cancer therapy using particle beams.
%
%    Copyright (C) 2014-2015 Particle Therapy Group Bergen
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
  names = {
      'histogram_uncertainty',
      'plateau_threshold',
      'linexp_alphas',
      'competition_params',
      'linplat_deltas'
    };
else
  names = varargin;
end

% Histogram uncertainties
% Estimated from the binning and significant figures in the DVH file data.
histogram_uncertainty = {
%   dose binning    Data point (y)
%   uncertainty     uncertainty
       0.05,            1e-6
  };

% Plateau Hall parameters per organ:
plateau_threshold = {
%                            Uncertainty
%    Organ name    Value    Low     High
    {'Stomach',      4.5,  -0.5,    0.5},
    {'Colon',        4.5,  -0.5,    0.5},
    {'Bone',         4.5,  -0.5,    0.5},
    {'Liver',        4.5,  -0.5,    0.5},
    {'Lungs',        4.5,  -0.5,    0.5},
    {'Bladder',      4.5,  -0.5,    0.5},
    {'Thyroid',      4.5,  -0.5,    0.5},
    {'Prostate',     4.5,  -0.5,    0.5},
    {'Breast',       4.5,  -0.5,    0.5},
  };

% Linear-exponential model alpha parameters per organ:
% Taken from table 1, Schneider et. al. 2005 paper:
%  "Estimation of radiation-induced cancer from three-dimensional dose
%   distributions: Concept of organ equivalent dose"
%  International journal of radiation oncology, biology, physics,
%  1 April 2005 (volume 61 issue 5 Pages 1510-1515
%  DOI: 10.1016/j.ijrobp.2004.12.040)
linexp_alphas = {
%                              Uncertainty
%    Organ name    Value     Low       High
    {'Stomach',    0.149,   -0.013,   0.014},
    {'Colon',      0.240,   -0.061,   0.075},
    {'Bone',       0.033,   -0.027,   0.043},
    {'Liver',      0.487,   -0.252,   0.385},
    {'Lungs',      0.129,   -0.012,   0.016},
    {'Bladder',    1.592,   -0.356,   0.434},
    {'Thyroid',    0.033,   -0.013,   0.014},
    {'Prostate',   0.804,   -0.496,   0.622},
    {'Breast',     0.041,   -0.041,   0.041},
  };
%breast parameters taken from Schneider 2011, no CI available

% Competition model alpha and alpha/beta ratio parameters per organ:
competition_params = {
%                            Uncertainty              Uncertainty   alpha/beta   Uncertainty
%    Organ name   alpha1    Low     High   alpha2    Low     High     ratio      Low    High
    {'Colon',     0.017,  -0.001,  0.001,   0.25,   -0.01,   0.01,     5.4,     -0.1,   0.1},
    {'Lungs',     0.017,  -0.001,  0.001,   0.25,   -0.01,   0.01,     4.5,     -0.1,   0.1},
    {'Bladder',   0.006,  -0.001,  0.001,   0.25,   -0.01,   0.01,     7.5,     -0.1,   0.1},
  };

% Linear-plateau delta parameters per organ:
linplat_deltas = {
%                              Uncertainty
%    Organ name    Value     Low       High
    {'Bladder',     0.1,    -0.01,     0.01},
  };


result = struct;
for n = 1:length(names)
  name = names{n};
  if ~ ischar(name)
    error('Expected a string for input parameter %d.', n);
  end
  switch name
    case 'histogram_uncertainty'
      result.dose_binning_uncertainty = histogram_uncertainty{1};
      result.volume_ratio_uncertainty = histogram_uncertainty{2};
    case 'plateau_threshold'
      for k = 1:length(plateau_threshold)
        p = plateau_threshold{k};
        result.(p{1}).plateau_threshold = struct('value', p{2}, 'range_low', p{2}+p{3}, 'range_high', p{2}+p{4});
      end
    case 'linexp_alphas'
      for k = 1:length(linexp_alphas)
        p = linexp_alphas{k};
        result.(p{1}).linexp_alpha = struct('value', p{2}, 'range_low', p{2}+p{3}, 'range_high', p{2}+p{4});
      end
    case 'competition_params'
      for k = 1:length(competition_params)
        p = competition_params{k};
        result.(p{1}).competition_alpha1 = struct('value', p{2}, 'range_low', p{2}+p{3}, 'range_high', p{2}+p{4});
        result.(p{1}).competition_alpha2 = struct('value', p{5}, 'range_low', p{5}+p{6}, 'range_high', p{5}+p{7});
        result.(p{1}).competition_alpha_beta_ratio = struct('value', p{8}, 'range_low', p{8}+p{9}, 'range_high', p{8}+p{10});
      end
    case 'linplat_deltas'
      for k = 1:length(linplat_deltas)
        p = linplat_deltas{k};
        result.(p{1}).linplat_delta = struct('value', p{2}, 'range_low', p{2}+p{3}, 'range_high', p{2}+p{4});
      end
    otherwise
      error('Unknown parameter set name %s.', name);
  end
end
