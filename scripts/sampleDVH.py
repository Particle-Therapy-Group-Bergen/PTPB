#!/usr/bin/env python

# This script provides a tool to sample DVHs and bootstrap them across patients.

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#    Particle Therapy Project Bergen (PTPB) - tools and models for research in
#    cancer therapy using particle beams.
#
#    Copyright (C) 2015 Particle Therapy Group Bergen
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Author: Artur Szostak <artursz@iafrica.com>

from __future__ import print_function
import os
import sys
import argparse
import subprocess
import textwrap
from numpy import arange, linspace
from processPatients import RunError, UncertaintyModel, Delta, DoubleDelta, \
                            Box, Box95, Triangle, Triangle95, Triangle95mode, \
                            Gaus, Gaus95, LogNorm, LogNorm95


class ConfigParams(object):
    """
    Object to represent configuration parameters loaded from the external
    configuration file.
    """

    def __init__(self, volume_bins, interpolation_method,
                 dose_binning_uncertainty, volume_ratio_uncertainty,
                 bootstrap_max_samples, bootstrap_sample_mode,
                 organ_name_map, organs):
        self.interpolation_method = interpolation_method
        self.dose_binning_uncertainty = dose_binning_uncertainty
        self.volume_ratio_uncertainty = volume_ratio_uncertainty
        self.bootstrap_max_samples = bootstrap_max_samples
        self.bootstrap_sample_mode = bootstrap_sample_mode
        self.organ_name_map = organ_name_map
        self.organs = organs
        self.volume_bins = volume_bins

    def __repr__(self):
        return "volume_bins = {0}\n" \
               "interpolation_method = {1}\n" \
               "dose_binning_uncertainty = {2}\n" \
               "volume_ratio_uncertainty = {3}\n" \
               "bootstrap_max_samples = {4}\n" \
               "bootstrap_sample_mode = {5}\n" \
               "organ_name_map = {6}\n" \
               "organs = {7}".format(
                   self.volume_bins, self.interpolation_method,
                   self.dose_binning_uncertainty, self.volume_ratio_uncertainty,
                   self.bootstrap_max_samples, self.bootstrap_sample_mode,
                   self.organ_name_map, self.organs)

    def generate_matlab_params(self):
        """generate_matlab_params() -> string

        Generates a string containing Matlab code to prepare the parameters
        needed for the sampleDVH.m function.
        """
        script = ""
        if self.volume_bins:
            script += "volumebins = {0};\n".format(self.volume_bins)
        # Prepare the organ_name_map variable.
        func = lambda x: "'{0}', '{1}'".format(x, self.organ_name_map[x])
        liststr = "; ".join(map(func, self.organ_name_map))
        organmapstr = "{" + liststr + "}"
        # Prepare the organ list.
        liststr = ", ".join(map(lambda x: "'{0}'".format(x), self.organs))
        organliststr = "{" + liststr + "}"
        script += textwrap.dedent("""\
            organ_name_map = {0};
            organs = {1};
            params = struct;
            params.interpolation_method = '{2}';
            params.bootstrap_max_samples = {3};
            params.bootstrap_sample_mode = '{4}';
            """.format(organmapstr, organliststr, self.interpolation_method,
                       self.bootstrap_max_samples, self.bootstrap_sample_mode))
        if self.dose_binning_uncertainty:
            script += textwrap.dedent("""\
                params.dose_binning_uncertainty_model = '{0}';
                params.dose_binning_uncertainty = {1};
                """.format(self.dose_binning_uncertainty.get_name(),
                           self.dose_binning_uncertainty.get_range()))
        else:
            script += textwrap.dedent("""\
                params.dose_binning_uncertainty_model = 'box';
                params.dose_binning_uncertainty = getParameters('histogram_uncertainty').dose_binning_uncertainty;
                """)
        if self.volume_ratio_uncertainty:
            script += textwrap.dedent("""\
                params.volume_ratio_uncertainty_model = '{0}';
                params.volume_ratio_uncertainty = {1};
                """.format(self.volume_ratio_uncertainty.get_name(),
                           self.volume_ratio_uncertainty.get_range()))
        else:
            script += textwrap.dedent("""\
                params.volume_ratio_uncertainty_model = 'box';
                params.volume_ratio_uncertainty = getParameters('histogram_uncertainty').volume_ratio_uncertainty;
                """)
        return script


def prepare_argument_parser():
    """prepare_argument_parser() -> object

    Prepares the command line argument parser.
    """
    argparser = argparse.ArgumentParser(
        description = """Samples patient DVHs by performing Monte-Carlo sampling
            of the DVH data, boot-strapping the values across input files and
            creating mean samples.""")
    argparser.add_argument("filelist", metavar = "<file>", nargs = '+',
        help = """DVH input files to process. Must be in Matlab format (.mat)
            that can be converted with the convertDVHtoMatlabFile.sh script from
            text files.""")
    argparser.add_argument("-c", "--config", dest = "configfile",
        default = None, metavar = "<file>", action = "store",
        help = """Provides a configuration file for the sampling parameters.
            The file uses python syntax.""")
    argparser.add_argument("-p", "--print", dest = "printconfig",
        default = False, action = "store_true",
        help = """Prints the configuration loaded.""")
    argparser.add_argument("-o", "--organ", dest = "organlist",
        default = [], metavar = "<name>", action = "append",
        help = """The name of an organ for which to process the DVH.""")
    argparser.add_argument("-O", "--outfile", dest = "outputfile",
        default = "output.mat", metavar = "<file>", action = "store",
        help = """The name of the output file to write.""")
    argparser.add_argument("-N", "--nsamples", dest = "nsamples",
        default = 10, metavar = "<number>", action = "store",
        help = """The number of samples to generate when performing the
            Monte-Carlo sampling.""")
    argparser.add_argument("-B", "--bins", dest = "bins",
        default = None, metavar = "<range>", action = "store",
        help = """The volume-ratio binning points for interpolation. Must be
            in Matlab colon notation for ranges.""")
    return argparser


def load_config(filename, organlist):
    """load_config(string, list) -> object

    Load parameter configurations from the given filename. Also need to provide
    the default organ list.
    """
    variables = {'Delta': Delta, 'DoubleDelta': DoubleDelta,
                 'Box': Box, 'Box95': Box95, 'Triangle': Triangle,
                 'Triangle95': Triangle95, 'Triangle95mode': Triangle95mode,
                 'Gaus': Gaus, 'Gaus95': Gaus95, 'LogNorm': LogNorm,
                 'LogNorm95': LogNorm95, 'arange': arange, 'linspace': linspace}
    try:
        execfile(filename, variables)
        interpolation_method = variables['interpolation_method']
        if not isinstance(interpolation_method, str):
            raise RunError("'interpolation_method' in the configuration file"
                           " '{0}' must be a string.".format(filename))
        dose_binning_uncertainty = variables['dose_binning_uncertainty']
        if not isinstance(dose_binning_uncertainty, UncertaintyModel):
            msg = "'dose_binning_uncertainty' in the configuration file '{0}'" \
                  " must be an uncertainty model object.".format(filename)
            raise RunError(msg)
        volume_ratio_uncertainty = variables['volume_ratio_uncertainty']
        if not isinstance(volume_ratio_uncertainty, UncertaintyModel):
            msg = "'volume_ratio_uncertainty' in the configuration file '{0}'" \
                  " must be an uncertainty model object.".format(filename)
            raise RunError(msg)
        bootstrap_max_samples = variables['bootstrap_max_samples']
        if not isinstance(bootstrap_max_samples, int):
            raise RunError("'bootstrap_max_samples' in the configuration file"
                           " '{0}' must be an integer.".format(filename))
        bootstrap_sample_mode = variables['bootstrap_sample_mode']
        if not isinstance(bootstrap_sample_mode, str):
            raise RunError("'bootstrap_sample_mode' in the configuration file"
                           " '{0}' must be a string.".format(filename))
        # The following fields are made optional.
        if 'volume_bins' in variables:
            volume_bins = variables['volume_bins']
            T = type(arange(0,1))
            if not isinstance(volume_bins, list) \
               and not isinstance(volume_bins, T):
                msg = "'volume_bins' in the configuration file '{0}' must be" \
                      " a list or numpy array.".format(filename)
                raise RunError(msg)
        else:
            volume_bins = None;
        if 'organ_name_map' in variables:
            organ_name_map = variables['organ_name_map']
        else:
            organ_name_map = {}
        if not isinstance(organ_name_map, dict):
            raise RunError("'organ_name_map' in the configuration file '{0}'"
                           " must be a dictionary.".format(filename))
        if 'organs' in variables:
            organs = variables['organs']
            if not isinstance(organ, list):
                raise RunError("'organs' in the configuration file '{0}' must"
                               " be a list of strings.".format(filename))
        else:
            organs = []
    except IOError as e:
        raise RunError(str(e))
    except KeyError as e:
        raise RunError("Missing declaration of variable {0} in the"
                       " configuration file '{1}'.".format(e, filename))
    organs += organlist
    organs = list(set(organs))  # make unique
    return ConfigParams(volume_bins, interpolation_method,
                        dose_binning_uncertainty, volume_ratio_uncertainty,
                        bootstrap_max_samples, bootstrap_sample_mode,
                        organ_name_map, organs)


def run():
    """run() -> None

    Parse the command line arguments and then process the patient DVH files.
    """
    argparser = prepare_argument_parser()
    args = argparser.parse_args()
    # Load the configuration file and print it if so requested.
    if args.configfile:
        params = load_config(args.configfile, args.organlist)
        if args.bins:
            # Override the volume bins from the command line.
            params.volume_bins = args.bins
    else:
        organs = list(set(args.organlist))  # make unique
        params = ConfigParams(args.bins, 'pchip', None, None, 6435, 'adaptive',
                              {}, organs)
    if args.printconfig:
        print("Parameter configuration:")
        if args.configfile:
            print(params)
        else:
            print("(No config given. Using default internal parameters.)")
    # Prepare the Matlab script to execute, starting with the DVH file list and
    # output file name.
    script = ""
    fileliststr = ", ".join(map(lambda x: "'{0}'".format(x), args.filelist))
    script += "dvhfiles = {" + fileliststr + "};\n"
    script += "outfilename = '{0}';\n".format(args.outputfile)
    script += "nsamples = {0};\n".format(args.nsamples)
    # Setup additional input parameters.
    script += params.generate_matlab_params()
    # Add commands to load the DVHs and collect the set of all organs.
    script += textwrap.dedent("""\
        dvhs = {};
        organset = struct;
        for n = 1:length(dvhfiles)
            dvhs{n} = getDoseVolumeHistogram(dvhfiles{n}, organ_name_map, params);
            organnames = fieldnames(dvhs{n});
            for m = 1:length(organnames)
                organset.(organnames{m}) = 1;
            end
        end
        if length(organs) == 0
            organs = fieldnames(organset);
        end
        """)
    # Add the function calls to merge data if the output file already exists.
    if os.path.exists(args.outputfile):
        script += textwrap.dedent("""\
            try
                Samples = load(outfilename, 'Samples').Samples;
            catch
                error(sprintf('Missing Samples in "%s".', outfilename));
            end
            if ~ isstruct(Samples)
                error(sprintf('Samples in "%s" is not a structure.', outfilename));
            end
            if ~ exist('volumebins')
                volumebins = Samples.volumebins;
            end
            if length(volumebins) ~= length(Samples.volumebins)
                error('The volume bins vector in "%s" has a different size than the specified binning.', outfilename);
            end
            if find(volumebins ~= Samples.volumebins)
                error('The volume bins vector in "%s" has different elements than the specified binning.', outfilename);
            end
            """)
    else:
        script += textwrap.dedent("""\
            if ~ exist('volumebins')
                volumebins = 0:0.1:1;
            end
            Samples = struct('doses', struct(), 'volumebins', volumebins);
            """)
    # Add the function call to generate/merge samples and save to file.
    script += textwrap.dedent("""\
        for n = 1:length(organs)
            organ = organs{n};
            dvh_subset = {};
            k = 1;
            for m = 1:length(dvhs)
                if isfield(dvhs{m}, organ)
                    dvh_subset{k} = dvhs{m}.(organ);
                    k += 1;
                end
            end
            if length(dvh_subset) == 0
                warning(sprintf('No DVHs found for organ "%s".', organ));
                continue;
            end
            S = sampleDVH(dvh_subset, nsamples, volumebins, params);
            if isfield(Samples.doses, (organ))
                Samples.doses.(organ) = [Samples.doses.(organ); S];
            else
                Samples.doses.(organ) = S;
            end
        end
        save('-7', outfilename, 'Samples');
        """)
    # Invoke octave to execute the Matlab script.
    command_less_script = ['octave', '-q', '--path', '@@MFILE_PATH@@', '--eval']
    command = command_less_script + [script]
    if subprocess.call(command) != 0:
        cmdstr = " ".join(command_less_script + ["'...'"])
        filename = "failed_script.m"
        with open(filename, "w") as scriptfile:
            scriptfile.write(script)
        raise RunError("Execution of command failed: {0}\nWrote the failing"
                       " script that to file '{1}'.".format(cmdstr, filename))


if __name__ == '__main__':
    try:
        run()
    except RunError as e:
        print("ERROR: " + str(e), file = sys.stderr)
