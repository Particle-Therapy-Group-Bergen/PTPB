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

# This script generates dose-volume distribution plots with gnuplot.

# Authors: Artur Szostak <artursz@iafrica.com>

# Function to print the command usage help message. Takes zero parameters.
function PrintHelp {
cat <<EOF
Usage: plotDoseVolumeDistribution.sh [options] filename [filename filename ...]

Creates a plot of dose volume distributions for organs from given DVH files. The
DVH files must first be converted with the convertDVHtoMatlabFile.sh tool to
Matlab (.mat) format. This tool will produce an intermediate gnuplot (.gp)
script and appropriate text data (.dat) files used by the gnuplot script to make
the plot in extended Postscript (.eps) format. If an output file is specified
with a different extension than .eps then the ImageMagick convert utility is
used to try convert the image to the given format.

The available options are:
    --help | -h
        Prints this help message.

    --output=<name> | -o=<name>
        Specifies the output file name to use for the plot. If the extension is
        not .eps then, the file is converted with the ImageMagick 'convert'
        utility, assuming that utility is available.

    --force | -f
        If indicated then any existing files are overwritten without prompting
        the user.

    --organ=<name> | -g=<name>
        Specifies an organ for which to plot the dose volume distribution.
        Multiple values can be given. If no values are given then the script
        tries to process all organs found in the data file.

    --verbose | -v
        If specified, verbose information messages are printed as the files are
        processed. Specify multiple times on the command line to get more
        verbose output.

This is a shell script that can be run from outside octave.

Running example: plotDoseVolumeDistribution.sh --output=plot.eps DVHfile.mat

EOF
}

# Generates Matlab commands to execute on the fly.
function MakeMatlabScript {
cat << EOF
organlist = {$ORGANS};
filenames = {$INPUT_FILES};
organs = {};
for m = 1:length(filenames)
  filename = filenames{m};
  printf('%s => %s and %s\n', filename, '$DATA_FILE', '$GNUPLOT_SCRIPT');
  data = load(filename, 'DVH_data');
  for n = 1:length(data.DVH_data.structures)
    s = data.DVH_data.structures{n};

    % Check if the organ is in the required list, if such a list was given.
    if length(organlist) > 0
      notfound = 1;
      for m = 1:length(organlist)
        if strcmp(s.structName, organlist{m})
          notfound = 0;
          break;
        end
      end
      if notfound
        continue;
      end
    end

    s.('filename') = filename;
    organs{length(organs)+1} = s;
  end
end

datestring = strftime("%c", localtime(time));

if length(organs) > 0
  % Find the organ with the longest dose array.
  reforgan = organs{1};
  for n = 1:length(organs)
    if length(organs{n}.dose) > length(reforgan.dose)
      reforgan = organs{n};
    end
  end
  maxlen = length(reforgan.dose);

  % Check the data ranges and units match.
  for n = 1:length(organs)
    if maxlen ~= length(organs{n}.dose)
      dlen = length(organs{n}.dose);
      if sum( reforgan.dose(1:dlen) - organs{n}.dose(1:dlen) ) ~= 0
        error('The dose data points for "%s" in "%s" and "%s" in "%s" are not compatible.',
              reforgan.structName, reforgan.filename, organs{n}.structName, organs{n}.filename);
      end
      % Fill with zeros to fill out to maxlen.
      organs{n}.dose = reforgan.dose;
      organs{n}.ratioToTotalVolume = [organs{n}.ratioToTotalVolume, zeros(1, maxlen - dlen)];
    end
    if length(organs{n}.dose) ~= length(organs{n}.ratioToTotalVolume)
      error('The number of dose data points and volume ratio values for "%s" in "%s" are different or missing.',
            organs{n}.structName, organs{n}.filename);
    end
    if ~ strcmp(reforgan.doseUnit, organs{n}.doseUnit)
      error('The dose data point units for "%s" in "%s" and "%s" in "%s" are different.',
            reforgan.structName, reforgan.filename, organs{n}.structName, organs{n}.filename);
    end
    if ~ strcmp(reforgan.ratioToTotalVolumeUnit, organs{n}.ratioToTotalVolumeUnit)
      error('The ratio to volume data point units for "%s" in "%s" and "%s" in "%s" are different.',
            reforgan.structName, reforgan.filename, organs{n}.structName, organs{n}.filename);
    end
  end

  % Write the data table to file.
  [file, msg] = fopen('$DATA_FILE', 'w');
  if file == -1
    error('Could not open the file "%s" for writing: %s', '$DATA_FILE', msg);
  end
  fprintf(file, '# Generated by plotDoseVolumeDistribution.sh on %s\n', datestring);
  fprintf(file, '#%-24s', ' Filename:');
  for n = 1:length(organs)
    fprintf(file, '%25s', organs{n}.filename);
  end
  fprintf(file, '\n');
  fprintf(file, '#%24s', ' ');
  for n = 1:length(organs)
    fprintf(file, '%25s', organs{n}.structName);
  end
  fprintf(file, '\n');
  heading = sprintf('%s [%s]', 'dose', reforgan.doseUnit);
  fprintf(file, '#%24s', heading);
  for n = 1:length(organs)
    heading = sprintf('%s [%s]', 'volume ratio', organs{n}.ratioToTotalVolumeUnit);
    fprintf(file, '%25s', heading);
  end
  fprintf(file, '\n');
  for m = 1:length(reforgan.dose)
    fprintf(file, '%25.16e', reforgan.dose(m));
    for n = 1:length(organs)
      fprintf(file, '%25.16e', organs{n}.ratioToTotalVolume(m));
    end
    fprintf(file, '\n');
  end
  if fclose(file) ~= 0
    error('There was a problem closing the file "%s".', '$DATA_FILE');
  end

  % Write the gnuplot script to file.
  [file, msg] = fopen('$GNUPLOT_SCRIPT', 'w');
  if file == -1
    error('Could not open the file "%s" for writing: %s', '$GNUPLOT_SCRIPT', msg);
  end
  fprintf(file, '# Generated by plotDoseVolumeDistribution.sh on %s\n', datestring);
  fprintf(file, 'set terminal postscript eps color enhanced dashed size 10cm,7.33cm\n');
  fprintf(file, 'set output "%s"\n', '$EPS_FILE');
  fprintf(file, 'set style data lines\n');
  fprintf(file, 'set title "Dose Volume Distribution"\n');
  fprintf(file, 'set xlabel "Dose (%s)"\n', reforgan.doseUnit);
  fprintf(file, 'set ylabel "Volume fraction (%s)"\n', reforgan.ratioToTotalVolumeUnit);
  fprintf(file, 'set xrange [%g:%g]\n', min(reforgan.dose), max(reforgan.dose));
  fprintf(file, 'set yrange [%g:%g]\n', min(reforgan.ratioToTotalVolume), max(reforgan.ratioToTotalVolume) * 1.1);
  fprintf(file, '#set key at 20, 80 samplen 0.8 spacing 1.25 font ",16"\n');
  fprintf(file, '#set logscale y\n');
  fprintf(file, '#set format y "10^{%%L}"\n');
  fprintf(file, '#set style line 1\n');
  fprintf(file, '#set size 0.9,0.9\n');
  for n = 1:length(organs)
    if n == 1
      fprintf(file, 'plot "%s"  using 1:%d  title "%s - %s"  lw 1',
              '$DATA_FILE', n+1, organs{n}.structName, organs{n}.filename);
    else
      fprintf(file, ', %s\n     "%s"  using 1:%d  title "%s - %s"  lw 1',
              '\\', '$DATA_FILE', n+1, organs{n}.structName, organs{n}.filename);
    end
  end
  fprintf(file, '\n');
  if fclose(file) ~= 0
    error('There was a problem closing the file "%s".', '$GNUPLOT_SCRIPT');
  end
end
EOF
}


# The path of the Matlab .m files. It is setup by the "Makefile install" target.
MFILE_PATH="@@MFILE_PATH@@"

# Set some initial defaults.
OUTPUT_FILENAME=""
FORCE_OVERWRITE=no
VERBOSE_COUNT=0
FILE_COUNT=0
ORGANS=""
INPUT_FILES=""

# Parse the command line parameters.
for ARG ; do
  case $ARG in
    --help|-h)
      PrintHelp
      exit
      ;;
    --output=*|-o=*)
      OUTPUT_FILENAME=`echo "$ARG" | sed 's|-.*=||'`;
      if test -z "$OUTPUT_FILENAME" ; then
        echo "ERROR: No output file name given for option '$ARG'." 1>&2
        exit 2
      fi
      ;;
    --force|-f)
      FORCE_OVERWRITE=yes
      ;;
    --organ=*|-g=*)
      ORGAN=`echo "$ARG" | sed 's|-.*=||'`
      if test -z "$ORGAN" ; then
        echo "ERROR: No organ name given for option '$ARG'." 1>&2
        exit 2
      fi
      if test -z "$ORGANS" ; then
        ORGANS="'$ORGAN'"
      else
        ORGANS="$ORGANS,'$ORGAN'"
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
      if test -z "$INPUT_FILES" ; then
        INPUT_FILES="'$ARG'"
      else
        INPUT_FILES="$INPUT_FILES,'$ARG'"
      fi
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


# Prepare the file name for the data table, gnuplot script and EPS output file name.
if test -n "$OUTPUT_FILENAME" ; then
  BASE_NAME=`dirname $OUTPUT_FILENAME`/`basename $OUTPUT_FILENAME | sed 's|\.[^\.]*$||'`
  DATA_FILE="${BASE_NAME}.dat"
  GNUPLOT_SCRIPT="${BASE_NAME}.gp"
  EPS_FILE="${BASE_NAME}.eps"
else
  DATA_FILE="plot.dat"
  GNUPLOT_SCRIPT="plot.gp"
  EPS_FILE="plot.eps"
fi

# Check if the gnuplot script or data tables already exists and skip it if either,
# of them do, unless the force option was used.
MAKE_SCRIPT_AND_TABLE=yes
if test -f "$DATA_FILE" ; then
  if test "$FORCE_OVERWRITE" = no ; then
    echo "WARNING: File '$DATA_FILE' already exists so skipping creation." 1>&2
    MAKE_SCRIPT_AND_TABLE=no
  fi
fi
if test -f "$GNUPLOT_SCRIPT" ; then
  if test "$FORCE_OVERWRITE" = no ; then
    echo "WARNING: File '$GNUPLOT_SCRIPT' already exists so skipping creation." 1>&2
    MAKE_SCRIPT_AND_TABLE=no
  fi
fi

if test "$MAKE_SCRIPT_AND_TABLE" = "yes" ; then
  # Generate the data table for the gnuplot script.
  if test "$VERBOSE_COUNT" -gt 1 ; then
    echo "============================== Matlab script ================================="
    MakeMatlabScript
    echo "=============================================================================="
  fi
  SCRIPT=`MakeMatlabScript`
  if test "$VERBOSE_COUNT" -gt 0 ; then
    octave $VERBOSE_OCTAVE --path "$MFILE_PATH" --eval "$SCRIPT" || exit $?
  else
    octave $VERBOSE_OCTAVE --path "$MFILE_PATH" --eval "$SCRIPT" > /dev/null || exit $?
  fi
fi

if test "$VERBOSE_COUNT" -gt 0 ; then
  echo "$GNUPLOT_SCRIPT and $DATA_FILE => $EPS_FILE"
fi

# Check if the output EPS file already exists and skip it if it does, unless the
# force option was used.
MAKE_EPS_FILE=yes
if test -f "$EPS_FILE" ; then
  if test "$FORCE_OVERWRITE" = no ; then
    echo "WARNING: File '$EPS_FILE' already exists so skipping creation." 1>&2
    MAKE_EPS_FILE=no
  fi
fi

if test "$MAKE_EPS_FILE" = "yes" ; then
  gnuplot $GNUPLOT_SCRIPT
fi

# Check if we can stop at the eps because we already got the file with gnuplot.
if test "`echo "$OUTPUT_FILENAME" | sed 's|^.*\(\.eps\)$|\1|'`" = ".eps" ; then
  exit
fi

# If we expect a different final output file format then use ImageMagick to
# convert the file. Note, we use a high dots per inch density to get good quality
# for raster formats like PNG or JPEG.
if test -n "$OUTPUT_FILENAME" ; then
  if test "$VERBOSE_COUNT" -gt 0 ; then
    echo "$EPS_FILE => $OUTPUT_FILENAME"
  fi
  if test -f "$OUTPUT_FILENAME" ; then
    if test "$FORCE_OVERWRITE" = no ; then
      echo "WARNING: File '$OUTPUT_FILENAME' already exists so skipping creation." 1>&2
      exit
    fi
  fi
  convert -density 600 "$EPS_FILE" "$OUTPUT_FILENAME"
fi
