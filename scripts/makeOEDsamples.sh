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

# This script is used to generate Organ Equivalent Dose (OED) calculation
# samples, by randomly sampling the undertainty distributions of the input model
# parameters and dose volume histogram data points.
# This allows estimation of propagated uncertainties using a monte carlo method.

# Authors: Artur Szostak <artursz@iafrica.com>

# Function to print the command usage help message. Takes zero parameters.
function PrintHelp {
cat <<EOF
Usage: makeOEDsamples.sh [options] filename [filename filename ...]

Calculates a number of samples of the Organ Equivalent Dose (OED) for dose
volume histograms found in a converted (.mat) DVH files. The files can be
converted with the convertDVHtoMatlabFile.sh script.
The sampling method is monte carlo sampling of the uncertainty distributions of
the model parameters and dose volume histogram data points. i.e. we jitter the
data points and parameters with a certain error distribution and then calculate
the OED each time. This allows one to perform uncertainty propagation through
the calculation using a monte carlo method.

The available options are:
    --help | -h
        Prints this help message

    --output=<file> | -o=<file>
        Specifies the full output file name to use.

    --suffix=<ext> | -s=<ext>
        Specifies the filename suffix and extension to use for the output file
        names. By default '_OED_results.mat' is used if this option is not
        given.

    --force | -f
        Will overwrite output files if they exit without prompting the user.

    --organ=<name> | -g=<name>
        Specifies an organ to process in the data file. Multiple values can be
        given. If no values are given then the script tries to process all
        organs found in the data file.

    --map=<key>:<value> | -a=<key>:<value>
        Indicates a mapping from a organ name found in the DVH file (the key) to
        a new standard name (value). The mapped name is used in the output and
        to search for response model parameters.

    --model=<name> | -m=<name>
        Gives the name of a response model to use for the OED calculation. More
        than one model can be given. If not models are explicitly passed on the
        command line then the script tries to compute all possible models.
        The list of valid models includes: 'LNT', 'PlateauHall', 'LinExp' or
        'Competition'.

    --samples=<number> | -n=<number>
        Indicates the number of samples to calculate.

    --dose-volume-uncertainty-model=<name> | -u=<name>
        Specifies the uncertainty distribution model to use for the dose volume
        histogram data points. This must be one of: 'box' or 'triangle'

    --parameter-uncertainty-model=<name> | -p=<name>
        Specifies the uncertainty distribution model to use for the response
        model parameters. This must be one of: 'box', 'triangle' or 'gaussian'.

    --integration-method=<name> | -i=<name>
        Indicates the integration method to use. Can be one of the following
        values: 'quad', 'quadv', 'quadl', quadgk' or 'trapz'.
        Refer to the documentation of the Matlab OED function for more details
        about the different methods. (Note some integration methods might not be
        supported with older versions of Matlab/Octave)

    --integration-tolerance=<number> | -t=<number>
        Specifies the tolerance threshold to use for the given integration
        method.

    --interpolation-method=<name> | -j=<name>
        Indicates the interpolation method to use. Can be one of the following
        values: 'nearest', 'linear', 'pchip', 'cubic' or 'spline'.
        Refer to the documentation of the Matlab interp1 function for more
        details about the different methods. (Note some interpolation methods
        might not be supported with older versions of Matlab/Octave)

    --dose-fractions=<number> |-k=<number>
        Gives the number of dose fractions to use in the competition model.

    --verbose | -v
        If specified, verbose information messages are printed as the files are
        converted. Specify multiple times on the command line to get more
        verbose output.

This is a shell script that can be run from outside octave.

Running example: makeOEDsamples.sh -v DVHinput.mat

EOF
}

# Generates Matlab commands to execute on the fly.
function MakeScript {
cat << EOF
options = struct;
organmap = {$ORGAN_MAP};
if length(organmap) > 0
  options.organ_name_map = organmap;
end
options.print_progress = $PRINT_PROGRESS;
options.debug_function = $DEBUG_FUNCTION;
options.struct_output = $STRUCT_OUTPUT;
options.dosevolume_uncertainty_model = '$DOSEVOLUME_UNCERTAINTY_MODEL';
options.parameter_uncertainty_model = '$PARAMETER_UNCERTAINTY_MODEL';
options.integration_method = '$INTEGRATION_METHOD';
options.integration_tolerance = $INTEGRATION_TOLERANCE;
options.interpolation_method = '$INTERPOLATION_METHOD';
options.dose_fractions = $DOSE_FRACTIONS;
OED_samples = sampleOED($NSAMPLES, '$INFILE', {$ORGANS}, {$MODELS}, options);
save('-7', '$OUTFILE', 'OED_samples');
EOF
}

# Generates Matlab commands to merge files.
function MakeMergeScript {
cat << EOF
try
  indata = load('$INFILE', 'OED_samples');
  outdata = load('$OUTFILE', 'OED_samples');
catch
  error("%s. \`OED_samples' structure not found in file.", lasterror.message);
end
OED_samples = outdata.OED_samples;
organs = fieldnames(indata.OED_samples);
for n = 1:length(organs)
  organname = organs{n}
  models = fieldnames(indata.OED_samples.(organname));
  for m = 1:length(models)
    modelname = models{m}
    x = indata.OED_samples.(organname).(modelname);
    OED_samples.(organname).(modelname) = [OED_samples.(organname).(modelname); x];
  end
end
save('-7', '$OUTFILE', 'OED_samples');
EOF
}


# The path of the Matlab .m files. It is setup by the "Makefile install" target.
MFILE_PATH="@@MFILE_PATHR@@"

# Set some initial defaults.
OUTPUT_FILENAME=""
OUTPUT_SUFFIX="_OED_results.mat"
FORCE_OVERWRITE=no
VERBOSE_COUNT=0
FILE_COUNT=0
NSAMPLES=10
ORGANS=""
ORGAN_MAP=""
MODELS=""
PRINT_PROGRESS=0
DEBUG_FUNCTION=0
DEBUG_FLAG=0
STRUCT_OUTPUT=1
DOSEVOLUME_UNCERTAINTY_MODEL="box"
PARAMETER_UNCERTAINTY_MODEL="box"
INTEGRATION_METHOD="trapz"
INTEGRATION_TOLERANCE=1e-4
INTERPOLATION_METHOD="linear"
DOSE_FRACTIONS=1
MERGE_RESULTS=no

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
    --suffix=*|-s=*)
      OUTPUT_SUFFIX=`echo "$ARG" | sed 's|-.*=||'`
      if test -z "$OUTPUT_SUFFIX" ; then
        echo "ERROR: No suffix given for option '$ARG'." 1>&2
        exit 2
      fi
      ;;
    --force|-f)
      FORCE_OVERWRITE=yes
      ;;
    --verbose|-v)
      let VERBOSE_COUNT++
      ;;
    --debug=*|-d=*)
      DEBUG_FLAG=`echo "$ARG" | sed 's|-.*=||'`
      ;;
    --debug|-d)
      DEBUG_FLAG=1
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
    --map=*|-a=*)
      KEYVAL=`echo "$ARG" | sed 's|-.*=||'`
      KEY=`echo "$KEYVAL" | sed 's|:.*$||'`
      VALUE=`echo "$KEYVAL" | sed 's|^.*:||'`
      if test -z "$KEY" ; then
        echo "ERROR: No key given for mapping option '$ARG'." 1>&2
        exit 2
      fi
      if test -z "$VALUE" ; then
        echo "ERROR: No value given for mapping option '$ARG'." 1>&2
        exit 2
      fi
      if test -z "$ORGAN_MAP" ; then
        ORGAN_MAP="'$KEY','$VALUE'"
      else
        ORGAN_MAP="$ORGAN_MAP;'$KEY','$VALUE'"
      fi
      ;;
    --model=*|-m=*)
      MODEL=`echo "$ARG" | sed 's|-.*=||'`
      if test -z "$MODEL" ; then
        echo "ERROR: No model name given for option '$ARG'." 1>&2
        exit 2
      fi
      if test -z "$MODELS" ; then
        MODELS="'$MODEL'"
      else
        MODELS="$MODELS,'$MODEL'"
      fi
      ;;
    --samples=*|-n=*)
      NSAMPLES=`echo "$ARG" | sed 's|-.*=||'`
      ;;
    --dose-volume-uncertainty-model=*|-u=*)
      DOSEVOLUME_UNCERTAINTY_MODEL=`echo "$ARG" | sed 's|-.*=||'`
      ;;
    --parameter-uncertainty-model=*|-p=*)
      PARAMETER_UNCERTAINTY_MODEL=`echo "$ARG" | sed 's|-.*=||'`
      ;;
    --integration-method=*|-i=*)
      INTEGRATION_METHOD=`echo "$ARG" | sed 's|-.*=||'`
      ;;
    --integration-tolerance=*|-t=*)
      INTEGRATION_TOLERANCE=`echo "$ARG" | sed 's|-.*=||'`
      ;;
    --interpolation-method=*|-j=*)
      INTERPOLATION_METHOD=`echo "$ARG" | sed 's|-.*=||'`
      ;;
    --dose-fractions=*|-k=*)
      DOSE_FRACTIONS=`echo "$ARG" | sed 's|-.*=||'`
      ;;
    --merge|-e)
      MERGE_RESULTS=yes
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

if test "$MERGE_RESULTS" = "yes" ; then
  if test -z "$OUTPUT_FILENAME" ; then
    OUTPUT_FILENAME="merged_OED_results.mat"
  fi
  if test -f "$OUTPUT_FILENAME" ; then
    if test "$FORCE_OVERWRITE" = no ; then
      echo "WARNING: File '$OUTPUT_FILENAME' already exists so skipping creation." 1>&2
      exit
    fi
  fi
else
  if test "$FILE_COUNT" -gt 1 -a -n "$OUTPUT_FILENAME" ; then
    echo "ERROR: Cannot use --option | -o if more than one input file is given." 1>&2
    exit 2
  fi
fi

VERBOSE_OCTAVE="--silent"
if test "$VERBOSE_COUNT" -gt 4 ; then
  VERBOSE_OCTAVE="--verbose"
fi
if test "$VERBOSE_COUNT" -gt 0 ; then
  let PRINT_PROGRESS=$VERBOSE_COUNT-1
fi

if test "$DEBUG_FLAG" -gt 1 ; then
  DEBUG_FUNCTION=1
fi


for ARG ; do
  case $ARG in
    # Ignore option strings
    --help|-h|--output=*|-o=*|--force|-f|--verbose|-v|--organ=*|-g=*) ;;
    --map=*|-a=*|--debug=*|-d=*|--debug|-d|--model=*|-m=*|--samples=*|-n=*) ;;
    --dose-volume-uncertainty-model=*|-u=*|--parameter-uncertainty-model=*|-p=*) ;;
    --integration-method=*|-i=*|--integration-tolerance=*|-t=*) ;;
    --interpolation-method=*|-j=*|--dose-fractions=*|-k=*|--merge|-e) ;;

    *)
      # Assume this argument is an input filename, so process it.
      SCRIPT=""
      INFILE="$ARG"
      if test -n "$OUTPUT_FILENAME" ; then
        OUTFILE="$OUTPUT_FILENAME"
      else
        OUTFILE="`dirname $INFILE`/`basename $INFILE | sed 's|\.[^\.]*$||'`$OUTPUT_SUFFIX"
      fi
      if test "$VERBOSE_COUNT" -gt 0 ; then
        echo "$INFILE => $OUTFILE"
      fi

      if test "$MERGE_RESULTS" = "yes" ; then
        # Setup script to merge multiple files into the output.
        if test -f "$OUTFILE" ; then
          if test "$DEBUG_FLAG" -gt 0 ; then
            echo "============================== Matlab script ================================="
            MakeMergeScript
            echo "=============================================================================="
          fi
          SCRIPT=`MakeMergeScript`
        else
          # Just copy the first file over. The next time it will be merged.
          cp "$INFILE" "$OUTFILE"
        fi
      else
        # Setup script to calculate OED values per file.
        if test -f "$OUTFILE" ; then
          if test "$FORCE_OVERWRITE" = no ; then
            echo "WARNING: File '$OUTFILE' already exists so skipping creation." 1>&2
            continue
          fi
        fi
        if test "$DEBUG_FLAG" -gt 0 ; then
          echo "============================== Matlab script ================================="
          MakeScript
          echo "=============================================================================="
        fi
        SCRIPT=`MakeScript`
      fi

      # Actually execute the script at this point.
      if test -n "$SCRIPT" ; then
        if test "$VERBOSE_COUNT" -gt 1 -o "$DEBUG_FUNCTION" -gt 0 ; then
          octave $VERBOSE_OCTAVE --path "$MFILE_PATH" --eval "$SCRIPT"
        else
          octave $VERBOSE_OCTAVE --path "$MFILE_PATH" --eval "$SCRIPT" > /dev/null
        fi
      fi
      ;;
  esac
done
