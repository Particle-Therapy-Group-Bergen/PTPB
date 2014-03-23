#!/bin/bash

###############################################################################
#
#    Particle Therapy Project Bergen (PTPB) - tools and models for research in
#    cancer therapy using particle beams.
#
#    Copyright (C) 2014 Particle Therapy Group Bergen
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
###############################################################################

# This script prints some statistics for sampled OED values found in a file.

# Authors: Artur Szostak <artursz@iafrica.com>

# Function to print the command usage help message. Takes zero parameters.
function PrintHelp {
cat <<EOF
Usage: printOEDstats.sh [options] filename [filename filename ...]

Prints a number of statistics for the Organ Equivalent Dose (OED) dose values
stored in the given files.

The available options are:
    --help | -h
        Prints this help message.

    --no-header | -n
        Indicates that textual headers should not be generated.

    --percent | -p
        Indicates that the standard deviation, min and max columns should be
        shown as a percentage from the mean value. The min and max columns will
        be the minimum deviation and maximum deviation from the mean in this
        case.

    --print-filename | -f
        If indicated then the file name is also printed in a column.

    --range=<matlab-range> | -r=<matlab-range>
        Used to indicate a range of values to process rather than all values.
        The <matlab-range> should be a Matlab vector syntax.

    --verbose | -v
        If specified, verbose information messages are printed as the files are
        processed. Specify multiple times on the command line to get more
        verbose output.

This is a shell script that can be run from outside octave.

Running example: printOEDstats.sh OEDresults.mat

EOF
}

# Generates Matlab commands to execute on the fly.
function MakeScript {
cat << EOF
make_header = ~ $NO_HEADER;
print_percent = $AS_PERCENT;
print_filename = $PRINT_FILENAME;
data = load('$INFILE', 'OED_samples');
if make_header
  if print_percent
    if print_filename
      printf('%12s%12s%12s%12s%12s%20s%20s  %s\n', 'mean', 'stdev %', 'min % dev', 'max % dev', '#samples', 'organ', 'model', 'filename');
    else
      printf('%12s%12s%12s%12s%12s%20s%20s\n', 'mean', 'stdev %', 'min % dev', 'max % dev', '#samples', 'organ', 'model');
    end
  else
    if print_filename
      printf('%12s%12s%12s%12s%12s%20s%20s  %s\n', 'mean', 'stdev', 'min', 'max', '#samples', 'organ', 'model', 'filename');
    else
      printf('%12s%12s%12s%12s%12s%20s%20s\n', 'mean', 'stdev', 'min', 'max', '#samples', 'organ', 'model');
    end
  end
end
% Count the number of rows.
count = 0;
organs = fieldnames(data.OED_samples);
for n = 1:length(organs)
  organ = organs{n};
  models = fieldnames(data.OED_samples.(organ));
  count += length(models);
end
% Now print the rows.
organs = fieldnames(data.OED_samples);
for n = 1:length(organs)
  organ = organs{n};
  models = fieldnames(data.OED_samples.(organ));
  for m = 1:length(models)
    model = models{m};
    x = data.OED_samples.(organ).(model);
    x = x($RANGE);
    count--;
    mean_x = mean(x);
    if print_percent
      std_x = std(x) / mean_x * 100;
      min_x = (min(x) - mean_x) / mean_x * 100;
      max_x = (max(x) - mean_x) / mean_x * 100;
    else
      std_x = std(x);
      min_x = min(x);
      max_x = max(x);
    end
    if (count > 0)
      if print_filename
        printf('%12g%12g%12g%12g%12d%20s%20s  %s\n', mean_x, std_x, min_x, max_x, length(x), organ, model, '$INFILE');
      else
        printf('%12g%12g%12g%12g%12d%20s%20s\n', mean_x, std_x, min_x, max_x, length(x), organ, model);
      end
    else
      if print_filename
        printf('%12g%12g%12g%12g%12d%20s%20s  %s', mean_x, std_x, min_x, max_x, length(x), organ, model, '$INFILE');
      else
        printf('%12g%12g%12g%12g%12d%20s%20s', mean_x, std_x, min_x, max_x, length(x), organ, model);
      end
    end
  end 
end
EOF
}


# The path of the Matlab .m files. It is setup by the "Makefile install" target.
MFILE_PATH="@@MFILE_PATHR@@"

# Set some initial defaults.
VERBOSE_COUNT=0
NO_HEADER=0
AS_PERCENT=0
PRINT_FILENAME=0
FILE_COUNT=0
RANGE="1:length(x)"

# Parse the command line parameters.
for ARG ; do
  case $ARG in
    --help|-h)
      PrintHelp
      exit
      ;;
    --no-header|-n)
      NO_HEADER=1
      ;;
    --percent|-p)
      AS_PERCENT=1
      ;;
    --print-filename|-f)
      PRINT_FILENAME=1
      ;;
    --range=*|-r=*)
      RANGE=`echo "$ARG" | sed 's|-.*=||'`
      if test -z "$RANGE" ; then
        echo "ERROR: No range given for option '$ARG'." 1>&2
        exit 2
      fi
      ;;
    --verbose|-v)
      let VERBOSE_COUNT++
      ;;
    *)
      # Check if the file exits:
      if test ! -f "$ARG" ; then
        echo "ERROR: The file '$ARG' could not be found." 1>&2
        exit 1
      fi
      let FILE_COUNT++
      ;;
  esac
done


# Check that we got at least one input file.
if test "$FILE_COUNT" -eq 0 ; then
  echo "ERROR: No input files given." 1>&2
  PrintHelp
  exit 1
fi

VERBOSE_OCTAVE="--silent"
if test "$VERBOSE_COUNT" -gt 2 ; then
  VERBOSE_OCTAVE="--verbose"
fi


for ARG ; do
  case $ARG in
    # Ignore option strings
    --help|-h|--verbose|-v|--no-header|-n|--percent|-p|--print-filename|-f|--range=*|-r=*)
      ;;
    *)
      # Assume this argument is an input filename, so process it.
      SCRIPT=""
      INFILE="$ARG"
      if test "$VERBOSE_COUNT" -gt 0 ; then
        echo "Processing $INFILE"
      fi
      if test "$VERBOSE_COUNT" -gt 1 ; then
        echo "============================== Matlab script ================================="
        MakeScript
        echo "=============================================================================="
      fi
      SCRIPT=`MakeScript`
      # Actually execute the script at this point.
      octave $VERBOSE_OCTAVE --path "$MFILE_PATH" --eval "$SCRIPT"
      ;;
  esac
done
