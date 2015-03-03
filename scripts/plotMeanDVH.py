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
import re
import math
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
        default = "output.eps", metavar = "<file>", action = SetUnique,
        help = """The name of the output file for the plot. If gnuplot script
            creation was requested then two additional files will be produced
            with the same basename but extensions '.gp' and '.dat'.""")
    argparser.add_argument("-Q", "--quantiles", dest = "quantiles",
        default = "0:0.1:1", metavar = "<range>", action = SetUnique,
        help = """The quantiles binning vector. Must be in Matlab colon notation
            for ranges or vector declaration syntax.""")
    argparser.add_argument("-e", "--extra-curves", dest = "num_extra_curves",
        default = 0, action = SetUnique, type = int,
        help = """Indicates the number of extra linearly interpolated curves to
            draw between each pair of quantile curves. The default is zero.""")
    argparser.add_argument("-B", "--bins", dest = "bins",
        default = "0:1:100", metavar = "<range>", action = SetUnique,
        help = """The dose binning points for interpolation. Must be in Matlab
            colon notation for ranges or vector declaration syntax.""")
    argparser.add_argument("-I", "--interpolation", dest = "interpMethod",
        default = "pchip", metavar = "<method>", action = SetUnique,
        help = """The interpolation method to use. See Octave interp1() for
            available values. The default is 'pchip'.""")
    argparser.add_argument("-D", "--density", dest = "plotdensity",
        default = False, action = "store_true",
        help = """If given then the output plot will be a uncertainty density
            plot rather than quantile curves. By default quantile curves are
            drawn.""")
    argparser.add_argument("-g", "--gnuplot", dest = "make_gnuplot",
        default = False, action = "store_true",
        help = """Will produce a gnuplot script and data file to produce
            publication quality plots.""")
    argparser.add_argument("-w", "--white-background", dest = "makewhitebkg",
        default = False, action = "store_true",
        help = """Will produce a plot with a white background. This is only
            relevant for density (-D) plots using the gnuplot script generation
            (-g).""")
    return argparser


def make_matlab_script(args):
    """make_matlab_script(object) -> string

    Prepare the Matlab script to execute.
    """
    script = ""
    # Setup the DVH file list.
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
    if args.make_gnuplot:
        filename, ext = os.path.splitext(args.outputfile)
        script += "outfilename = '{0}';\n".format(filename);
    else:
        script += "outfilename = '{0}';\n".format(args.outputfile);
    script += "interpMethod = '{0}';\n".format(args.interpMethod);
    script += "dosebins = {0};\n".format(args.bins);
    script += "quantilebins = {0};\n".format(args.quantiles);
    script += "num_extra_curves = {0};\n".format(args.num_extra_curves);
    script += "plotdensity = {0};\n".format('1' if args.plotdensity else '0');
    script += "make_gnuplot = {0};\n".format('1' if args.make_gnuplot else '0');
    # Add commands to load the DVHs.
    script += textwrap.dedent("""\
        data = [];  % The data matrix for gnuplot data table.
        data_organs = {};  % The organ names added to the data table.
        data_index = [];  % The column index where the n'th organ's data starts.
        data_N = [];  % The number of data columns related to the n'th organ.
        data_i = 1;
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
                % Interpolate the curves linearly along the Y direction to get a
                % higher curve density.
                [nr, nc] = size(Y);
                newY = [];
                lambda = linspace(0, 1, 2 + num_extra_curves);
                lambda = repmat(lambda(2:length(lambda)-1), nr, 1);
                for k = 1:nc-1
                    y1 = repmat(Y(:,k), 1, size(lambda, 2));
                    y2 = repmat(Y(:,k+1), 1, size(lambda, 2));
                    interpolatedY = y1 .* (1 - lambda) + y2 .* lambda;
                    newY = [newY, Y(:,k), interpolatedY];
                end
                newY = [newY, Y(:,nc)];
                Y = newY;
                X = repmat(X(:,1), 1, size(Y, 2));
                if make_gnuplot
                    if length(data) == 0
                        data = X(:,1);
                    end
                    data_organs{data_i} = sprintf('%s:%s', inputfile, organ);
                    data_index(data_i) = size(data, 2);
                    data_N(data_i) = size(Y, 2);
                    data_i += 1;
                    data = [data Y];
                else
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
        end
        if make_gnuplot
            save('-ascii', sprintf('%s.dat', outfilename), 'data');
            printf('%g %g\\n', min(data(:,1)), max(data(:,1)));
            for n = 1:data_i-1
                printf('%d %d %s\\n', data_index(n), data_N(n), data_organs{n});
            end
        else
            if plotdensity
                title('Uncertainty density for mean dose volume histogram');
            else
                title('Quantile contours for mean dose volume histogram');
            end
            xlabel('Dose');
            ylabel('Volume ratio');
            print('-landscape', '-color', outfilename);
        end
        """)
    return script


class OrganInfo(object):
    """
    Class to contain information about organ data written by Matlab script.
    """
    def __init__(self, index, ncols, idname):
        self.index = index    # Start index into written data table.
        self.ncols = ncols    # The number of column of data in the data table.
        self.idname = idname  # The ID string identifying the organ.


def make_gnuplot_density_script(args, organinfos, xmin, xmax):
    """make_gnuplot_density_script(object, list, float, float) ->
                                                            (string, string)

    Prepares the gnuplot script for producing a mean dose volume histogram
    density plot. Teturns the tuple (name, ext), which is the name of the script
    file written and the extension of the produced file.
    """
    basename, ext = os.path.splitext(args.outputfile)
    filename = basename + '.gp'
    script = textwrap.dedent("""\
        # The following are alternative targets for testing:
        #set terminal postscript eps color enhanced dashed size 8cm,6cm
        #set output "@@BASENAME@@.eps"
        #set terminal pdfcairo enhanced dashed size 8cm,6cm
        #set output "@@BASENAME@@.pdf"
        #set terminal pngcairo enhanced dashed transparent truecolor font "times,75" size 1890,1418 linewidth 1

        # Note: "font times 100" gives the equivalent of 8pt font in openoffice once rescaled to 8cm x 6cm.
        # The size is chosen to give 600 dpi resolution for a 8cm by 6cm plot.
        set terminal png enhanced transparent truecolor font times 100 size 3780,2836 linewidth 3
        set output "@@BASENAME@@.png"

        set style data lines
        set title "Mean dose volume histogram."
        set xlabel "Dose [Gy(RBE)]"
        set ylabel "Volume fraction [%]"
        set xrange [@@XMIN@@:@@XMAX@@]
        set yrange [0:100]
        set border linewidth 3

        #set key at 20, 80 samplen 0.8 spacing 1.25 font ",16"
        #set logscale y
        #set format y "10^{%L}"
        #set style line 1
        #set size 0.9,0.9

        plot \\
        """.replace("@@BASENAME@@", basename).replace(
                    "@@XMIN@@", xmin).replace(
                    "@@XMAX@@", xmax))
    color = ["#0000FF", "#00FF00", "#FF0000", "#00FFFF", "#FFFF00", "#FF00FF"]
    colorn = 0
    plotlines = []
    for info in organinfos:
        maxk = int(math.ceil(info.ncols * 0.5))
        for k in xrange(0, maxk):
            title = "notitle" if k > 0 else 'title "{0}"'.format(info.idname)
            plotline = '  "{0}.dat" using 1:(${1}*100):(${2}*100) {3}' \
                ' linecolor rgb "{4}" with filledcurves fs transparent' \
                ' solid {5}'.format(basename, info.index + k + 1,
                                    info.index + info.ncols - k, title,
                                    color[colorn], 0.999 / maxk)
            plotlines.append(plotline)
        colorn = (colorn + 1) % len(color)
    script += ", \\\n".join(plotlines)
    with open(filename, 'w') as scriptfile:
        scriptfile.write(script)
    return (filename, '.png')


def make_gnuplot_line_script(args, organinfos, xmin, xmax):
    """make_gnuplot_line_script(object, list, float, float) -> (string, string)

    Prepares the gnuplot script for producing an EPS line plot of the dose
    volume histogram. Returns the tuple (name, ext), which is the name of the
    script file written and the extension of the produced file.
    """
    basename, ext = os.path.splitext(args.outputfile)
    filename = basename + '.gp'
    script = textwrap.dedent("""\
        #set terminal pdfcairo enhanced dashed size 8cm,6cm
        #set output "@@BASENAME@@.pdf"
        #set terminal pngcairo enhanced dashed transparent truecolor font "times,75" size 1890,1418 linewidth 1
        # Note: "font times 100" gives the equivalent of 8pt font in openoffice once rescaled to 8cm x 6cm.
        # The size is chosen to give 600 dpi resolution for a 8cm by 6cm plot.
        #set terminal png enhanced transparent truecolor font times 100 size 3780,2836 linewidth 3
        #set output "@@BASENAME@@.png"

        set terminal postscript eps color enhanced dashed size 8cm,6cm
        set output "@@BASENAME@@.eps"

        set style data lines
        set title "Mean dose volume histogram."
        set xlabel "Dose [Gy(RBE)]"
        set ylabel "Volume fraction [%]"
        set xrange [@@XMIN@@:@@XMAX@@]
        set yrange [0:100]
        set border linewidth 1

        #set key at 20, 80 samplen 0.8 spacing 1.25 font ",16"
        #set logscale y
        #set format y "10^{%L}"
        #set style line 1
        #set size 0.9,0.9

        plot \\
        """.replace("@@BASENAME@@", basename).replace(
                    "@@XMIN@@", xmin).replace(
                    "@@XMAX@@", xmax))
    color = ["#0000FF", "#00FF00", "#FF0000", "#00FFFF", "#FFFF00", "#FF00FF"]
    colorn = 0
    plotlines = []
    for info in organinfos:
        maxk = int(math.ceil(info.ncols * 0.5))
        for k in xrange(0, info.ncols):
            title = "notitle" if k > 0 else 'title "{0}"'.format(info.idname)
            plotline = '  "{0}.dat" using 1:(${1}*100) {2} lw 1 lt 1' \
                ' linecolor rgb "{3}"'.format(basename, info.index + k + 1,
                                              title, color[colorn])
            plotlines.append(plotline)
        colorn = (colorn + 1) % len(color)
    script += ", \\\n".join(plotlines)
    with open(filename, 'w') as scriptfile:
        scriptfile.write(script)
    return (filename, '.eps')


def make_gnuplot_script(args, organinfos, xmin, xmax):
    """make_gnuplot_script(object, list, float, float) -> (string, string)

    Produces a density plot if args.plotdensity is True and a line plot
    otherwise.
    """
    if args.plotdensity:
        return make_gnuplot_density_script(args, organinfos, xmin, xmax)
    else:
        return make_gnuplot_line_script(args, organinfos, xmin, xmax)


def run():
    """run() -> None

    Parse the command line arguments and then process the files with DVH
    samples.
    """
    argparser = prepare_argument_parser()
    args = argparser.parse_args()
    if args.num_extra_curves < 0:
        raise RunError("The number of extra curves must be a positive integer.")
    script = make_matlab_script(args)
    # Invoke octave to execute the Matlab script.
    command_less_script = ['octave', '-q', '--path', '@@MFILE_PATH@@', '--eval']
    command = command_less_script + [script]
    proc = subprocess.Popen(command, stdout = subprocess.PIPE)
    out, err = proc.communicate()
    if proc.returncode != 0:
        cmdstr = " ".join(command_less_script + ["'...'"])
        filename = "failed_script.m"
        with open(filename, "w") as scriptfile:
            scriptfile.write(script)
        raise RunError("Execution of command failed: {0}\nWrote the failing"
                       " script that to file '{1}'.".format(cmdstr, filename))
    if args.make_gnuplot:
        # Parse the size out the data matrix.
        pattern1 = re.compile(r'^(\S+)\s+(\S+)$')
        pattern2 = re.compile(r'^(\d+)\s+(\d+)\s+(.*)$')
        organinfos = []
        for n, line in enumerate(out.splitlines()):
            match = None
            if n == 0:
                match = pattern1.match(line)
                if match is not None:
                    xmin = match.group(1)
                    xmax = match.group(2)
            else:
                match = pattern2.match(line)
                if match is not None:
                    index = int(match.group(1))
                    ncols = int(match.group(2))
                    organID = match.group(3)
                    info = OrganInfo(index, ncols, organID)
                    organinfos.append(info)
            if match is None:
                raise RunError("The output received from the octave command"
                               " does not appear to be valid:\n" + out)
        # Prepare the gnuplot script and execute it.
        scriptname, created_file_ext = make_gnuplot_script(args, organinfos,
                                                           xmin, xmax)
        if subprocess.call(['gnuplot', scriptname]) != 0:
            raise RunError("Execution of gnuplot command failed.")
        basename, required_ext = os.path.splitext(args.outputfile)
        gpout_filename = basename + created_file_ext
        # Add white background if requested.
        if args.plotdensity and args.makewhitebkg:
            if subprocess.call(['convert', gpout_filename, '-background',
                                'white', '-flatten', gpout_filename]) != 0:
                raise RunError("Execution of convert command failed.")
        # Convert to the output file format if not just a PNG is required.
        if created_file_ext != required_ext:
            if subprocess.call(['convert', '-density', '600', gpout_filename,
                                args.outputfile]) != 0:
                raise RunError("Execution of convert command failed.")


if __name__ == '__main__':
    try:
        run()
    except RunError as e:
        print("ERROR: " + str(e), file = sys.stderr)
