#!/usr/bin/env python

# This script is used to calculate how a dose response function changes a DVH.
# It will load a file produced by the convertDVHtoMatlabFile.sh script, apply
# the configured dose responses and write out the modified DVH's to file.

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
import re
import argparse
import collections
import subprocess
import textwrap
from processPatients import RunError, SetUnique


def prepare_argument_parser():
    """prepare_argument_parser() -> object

    Prepares the command line argument parser.
    """
    argparser = argparse.ArgumentParser(
        description = """Applies dose response models to one or more dose volume
            histograms (DVHs) loaded from file. The modified DVHs are then
            written out to file again in the same overall format as the input
            files.""")
    argparser.add_argument("filelist", metavar = "<file>", nargs = '+',
        help = """One or more input DVH data files.""")
    argparser.add_argument("outpath", metavar = "<path>", nargs = 1,
        help = """The output file name or directory if multiple input files are
            given.""")
    argparser.add_argument("-p", "--print-organs", dest = "print_organs",
        default = False, action = "store_true",
        help = """If specified the list of available organs will be printed and
            no files produced.""")
    argparser.add_argument("-o", "--organ", dest = "organlist",
        default = [], metavar = "<name>", action = "append",
        help = """The name of an organ to process. Note that this must be the
            structure name exactly as found in the DVH file. If no organ names
            are given then all structures are processed.""")
    argparser.add_argument("-m", "--model", dest = "model",
        default = "LinExp", metavar = "<name>", action = SetUnique,
        help = """The name of a response model or function to apply. Note that
            this can be a Matlab anonymous function declaration.""")
    argparser.add_argument("-a", "--arg", dest = "func_params",
        default = [], metavar = "<value>", action = "append",
        help = """An argument to pass to the response model function. More than
            one argument can be specified. They are passed to the function in
            the order given on the command line. The syntax for the arguments
            must be Matlab compatible.""")
    argparser.add_argument("-D", "--dose-bins", dest = "dosebins",
        default = "0:1:100", metavar = "<range>", action = SetUnique,
        help = """The dose binning points for interpolation. Must be in Matlab
            colon notation for ranges or vector declaration syntax.""")
    argparser.add_argument("-V", "--volume-ratio-bins", dest = "volumebins",
        default = "0:0.001:1", metavar = "<range>", action = SetUnique,
        help = """The volume ratio binning points for interpolation. Must be in
            Matlab colon notation for ranges or vector declaration syntax.""")
    argparser.add_argument("-I", "--interpolation", dest = "interpMethod",
        default = "pchip", metavar = "<method>", action = SetUnique,
        help = """The interpolation method to use. See Octave interp1() for
            available values. The default is 'pchip'.""")
    return argparser


def make_matlab_script(args):
    """make_matlab_script(object) -> string

    Prepare the Matlab script to execute based on the parsed arguments.
    """
    script = ""
    # Setup the DVH file list.
    fileliststr = ", ".join(map(lambda x: "'{0}'".format(x), args.filelist))
    script += "inputfiles = {" + fileliststr + "};\n"
    # Add the output file path
    script += "outputpath = '{0}';\n".format(args.outpath[0])
    # Setup the organ name list if it is not empty.
    if len(args.organlist) > 0:
        organliststr = ", ".join(map(lambda x: "'{0}'".format(x),
                                     args.organlist))
        script += "organlist = {" + organliststr + "};\n"
    # Setup the model function.
    pattern = re.compile(r'^[A-Za-z0-9_]+$')
    if pattern.match(args.model):
        script += "responseFunc = '{0}';\n".format(args.model)
    else:
        script += "responseFunc = {0};\n".format(args.model)
    # Add any response function arguments as a cell array.
    argstr = ", ".join(map(lambda x: str(x), args.func_params))
    script += "responseArgs = {" + argstr + "};\n"
    # Add the print_organs flag.
    script += "print_organs = {0};\n".format('1' if args.print_organs else '0')
    # Add the dose and volume ratio binning variables.
    script += "dosebins = {0};\n".format(args.dosebins)
    script += "volumebins = {0};\n".format(args.volumebins)
    # Add the interpolation method to use.
    script += "interpMethod = '{0}';\n".format(args.interpMethod);
    # Now add the main part of the script.
    script += textwrap.dedent("""\
        for n = 1:length(inputfiles)
            inputfile = inputfiles{n};
            try
                old_DVH_data = load(inputfile, 'DVH_data').DVH_data;
            catch
                error(sprintf('Missing DVH_data in "%s".', inputfile));
            end
            if ~ isstruct(old_DVH_data)
                error(sprintf('DVH_data in "%s" is not a structure.', inputfile));
            end
            % Get list of all structures in the input file.
            allStructures = {};
            for m = 1:length(old_DVH_data.structures)
                allStructures{m} = old_DVH_data.structures{m}.structName;
            end
            if print_organs
                for m = 1:length(allStructures)
                    printf('%s\\n', allStructures{m});
                end
                continue
            end
            if exist('organlist')
                % Check that the organs passed on the command line exist in the file.
                for m = 1:length(organlist)
                    if ~ any(ismember(allStructures, organlist{m}))
                        error('The organ "%s" does not exist in "%s".',
                              organlist{m}, inputfile);
                    end
                end
                organsToProcess = organlist;
            else
                organsToProcess = allStructures;
            end
            % Create a name mapping for getDoseVolumeHistogram so that we do not
            % by mistake use organ names that cannot be structure field names.
            namemap = {};
            nametoindex = struct;  % Use this as a mapping back to the structure number.
            for m = 1:length(organsToProcess)
                newname = sprintf('organ%d', m);
                namemap{m,1} = organsToProcess{m};
                namemap{m,2} = newname;
                nametoindex.(newname) = m;
            end
            renamed_organs = fieldnames(nametoindex);
            % Prepare DVH data structure.
            DVH_data = struct('header', old_DVH_data.header, 'structures', {{}});
            structnum = 1;
            % Apply the response function to the DVH's and append these to the
            % new DVH_data structure.
            data = getDoseVolumeHistogram(inputfile, namemap, renamed_organs{:});
            for m = 1:length(renamed_organs)
                organ = renamed_organs{m};
                dvh = data.(organ);
                newdvh = applyResponse(dvh, volumebins, dosebins, interpMethod,
                                       responseFunc, responseArgs{:});
                index = nametoindex.(organ);
                dose = newdvh.datapoints(:,1);
                vol = newdvh.datapoints(:,2);
                DVH_data.structures{structnum} = old_DVH_data.structures{index};
                DVH_data.structures{structnum}.dose = dose';
                DVH_data.structures{structnum}.ratioToTotalVolume = vol';
                structnum += 1;
            end
            % Save to the output file
            if length(inputfiles) == 1
                save('-7', outputpath, 'DVH_data');
            else
                [pathstr, fname, ext] = fileparts(inputfile);
                filename = strcat(fname, ext);
                save('-7', sprintf('%s/%s', outputpath, filename), 'DVH_data');
            end
        end
        """)
    return script


def run():
    """run() -> None

    Parse the command line arguments and then process the DVH data files.
    """

    def _call_command(*args, **kwargs):
        try:
            return subprocess.call(*args, **kwargs)
        except OSError as err:
            msg = "Failed to run '{0}': {1}".format(args[0][0], err)
            raise RunError(msg)

    argparser = prepare_argument_parser()
    args = argparser.parse_args()
    # Prepare the output path if more than one input file was given.
    if len(args.filelist) > 1:
        if not os.path.exists(args.outpath[0]):
            os.makedirs(args.outpath[0])
        if not os.path.isdir(args.outpath[0]):
            raise RunError("The output path '{0}' is not a"
                           " directory".format(args.outpath[0]))
    script = make_matlab_script(args)
    # Invoke octave to execute the Matlab script.
    command_less_script = ['octave', '-q', '--path', '@@MFILE_PATH@@', '--eval']
    command = command_less_script + [script]
    if _call_command(command) != 0:
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
