#!/usr/bin/env python

# This script produces mean DVH plots from sampled data.

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
from processPatients import RunError


def prepare_argument_parser():
    """prepare_argument_parser() -> object

    Prepares the command line argument parser.
    """
    argparser = argparse.ArgumentParser(
        description = """Plots the mean DVH and uncertainty distribution from
            sampled data produced by the sampleDVH.py tool.""")
    argparser.add_argument("filelist", metavar = "<file>", nargs = '+',
        help = """One or more input sampled data files.""")
    argparser.add_argument("-o", "--organ", dest = "organlist",
        default = [], metavar = "<name>[:<file>]", action = "append",
        help = """The name of an organ for which to plot. A file name can be
            given after the colon to indicate that the organ selection should
            apply to a specific file.""")
    argparser.add_argument("-O", "--outfile", dest = "outputfile",
        default = "output.eps", metavar = "<file>", action = "store",
        help = """The name of the output file for the plot.""")
    argparser.add_argument("-Q", "--quantiles", dest = "quantiles",
        default = "0:0.1:1", metavar = "<range>", action = "store",
        help = """The quantiles binning vector. Must be in Matlab colon notation
            for ranges or vector declaration syntax.""")
    argparser.add_argument("-B", "--bins", dest = "bins",
        default = "0:1:100", metavar = "<range>", action = "store",
        help = """The dose binning points for interpolation. Must be in Matlab
            colon notation for ranges or vector declaration syntax.""")
    argparser.add_argument("-I", "--interpolation", dest = "interpMethod",
        default = "pchip", metavar = "<method>", action = "store",
        help = """The interpolation method to use. See Octave interp1() for
            available values. The default is 'pchip'.""")
    argparser.add_argument("-D", "--density", dest = "plotdensity",
        default = False, action = "store_true",
        help = """If given then the output plot will be a uncertainty density
            plot rather than quantile curves. By default quantile curves are
            drawn.""")
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
    script += "inputfiles = {" + fileliststr + "};\n"
    # Add the organ list.
    organlist = collections.defaultdict(list)
    for item in args.organlist:
        organ, seperator, filename = item.partition(':')
        if filename:
            organlist[organ].append(filename)
        else:
            organlist[organ].extend(args.filelist)
    script += "organlist = struct;\n"
    for key in organlist:
        fileliststr = ", ".join(map(lambda x: "'{0}'".format(x),
                                    organlist[key]))
        script += "organlist." + key + " = {" + fileliststr + "};\n"
    # Add other needed parameters.
    script += "outfilename = '{0}';\n".format(args.outputfile);
    script += "interpMethod = '{0}';\n".format(args.interpMethod);
    script += "dosebins = {0};\n".format(args.bins);
    script += "quantilebins = {0};\n".format(args.quantiles);
    script += "plotdensity = {0};\n".format('1' if args.plotdensity else '0');
    # Add commands to load the DVHs.
    script += textwrap.dedent("""\
        for n = 1:length(inputfiles)
            inputfile = inputfiles{n};
            try
                Samples = load(inputfile, 'Samples').Samples;
            catch
                error(sprintf('Missing Samples in "%s".', inputfile));
            end
            if ~ isstruct(Samples)
                error(sprintf('Samples in "%s" is not a structure.', inputfile));
            end
            if length(fieldnames(organlist)) > 0
                names = fieldnames(organlist);
                organs = {};
                k = 1;
                % Find all organ names that refer to the current inputfile.
                for m = 1:length(names)
                    [index, errmsg] = cellidx(organlist.(names{m}), inputfile);
                    if index > 0
                        organs{k} = names{m};
                        k += 1;
                    end
                end
            else
                organs = fieldnames(Samples.doses);
            end
            for m = 1:length(organs)
                organ = organs{m};
                [X, Y] = calcDVHquantiles(Samples.doses.(organ), Samples.volumebins,
                                          dosebins, quantilebins, interpMethod);
                if plotdensity
                    colormap hot;
                    [X, Y, Z] = contoursToMatrix(X(:,1), Y');
                    surf(X, Y, Z);
                    shading interp;
                    view(2);
                else
                    plot(X, Y);
                end
                hold on;
            end
        end
        print('-landscape', '-color', outfilename);
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
