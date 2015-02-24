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
import collections
import subprocess
import textwrap


class RunError(Exception):
    """Exception class for run time errors handled by this script."""
    pass


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
    argparser.add_argument("-o", "--organ", dest = "organlist",
        default = [], metavar = "<name>", action = "append",
        help = """A name of an organ to process.""")
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


def run():
    """run() -> None

    Parse the command line arguments and then process the patient DVH files.
    """
    argparser = prepare_argument_parser()
    args = argparser.parse_args()
    # Prepare the Matlab script to execute, starting with the DVH file list.
    script = ""
    fileliststr = ", ".join(map(lambda x: "'{0}'".format(x), args.filelist))
    script += "dvhfiles = {" + fileliststr + "};\n"
    # Setup additional input parameters.
    script += "outfilename = '{0}';\n".format(args.outputfile)
    script += "organ_name_map = {'Bladder_P', 'Bladder'};\n"
    script += "organ = 'Bladder';\n"
    script += "params = struct;\n"
    script += "params.dose_binning_uncertainty_model = 'box';\n"
    script += "params.volume_ratio_uncertainty_model = 'box';\n"
    script += "params.interpolation_method = 'pchip';\n"
    script += "params.bootstrap_max_samples = 6435;\n"
    script += "params.bootstrap_sample_mode = 'adaptive';\n"
    script += "nsamples = {0};\n".format(args.nsamples)
    if args.bins:
        script += "volumebins = {0};\n".format(args.bins)
    # Add commands to load the DVHs.
    script += textwrap.dedent("""\
        dvhs = {};
        for n = 1:length(dvhfiles)
            dvhs{n} = getDoseVolumeHistogram(dvhfiles{n}, organ_name_map, organ).(organ);
        end
        """)
    # Add the function calls to merge data if the output file already exists.
    if os.path.exists(args.outputfile):
        script += textwrap.dedent("""\
            Samples = load(outfilename, 'Samples').Samples;
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
            Samples = struct('doses', [], 'volumebins', volumebins);
            """)
    script += textwrap.dedent("""\
        S = sampleDVH(dvhs, nsamples, volumebins, params);
        Samples.doses = [Samples.doses; S];
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
