#!/bin/bash

# Function to print the command usage help message. Takes zero parameters.
function PrintHelp {
cat <<EOF
Usage: convertDVHtoMatlabFile.sh [options] filename [filename filename ...]

Converts one or more DVH input files into Octave .mat files. These can be
loaded much more quickly by Octave than parsing the raw DVH files with a
Matlab script.
Each input file will be converted and written to a file with the same basename
as the input file, but with the extension changed to .mat by default. For
example, myfile.txt will converted and written to myfile.mat. 

The available options are:
    --help | -h
        Prints this help message

    --output=<file> | -o=<file>
        Specifies the full output file name to use. This option is only valid if
        just one input filename is given.

    --suffix=<ext> | -s=<ext>
        Specifies the filename extension (without the dot) to use for the output
        file names. By default 'mat' (.mat) is used if this option is not given.

    --force | -f
        Will overwrite output files if they exit without prompting the user.

    --verbose | -v
        If specified verbose information messages are printed as the files are
        converted. Specify multiple times on the command line to get more
        verbose output.
EOF
}


# Get the path of the readDVHfile.m script. This should be in the same location
# as this script. We use the 'which' command to locate the path of this script,
# and $0 will contain its name.

MFILE_PATH=`which $0`
MFILE_PATH=`dirname $MFILE_PATH`

# Run through the list of arguments once to parse the options.
OUTPUT_FILENAME=""
OUTPUT_SUFFIX=""
FORCE_OVERWRITE=no
VERBOSE_COUNT=0
FILE_COUNT=0
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

if test "$FILE_COUNT" -eq 0 ; then
    exit
fi

if test "$FILE_COUNT" -gt 1 -a -n "$OUTPUT_FILENAME" ; then
    echo "ERROR: Cannot use --option | -o if more than one input file is given." 1>&2
    exit 2
fi

# Assign some default values if not set by the command line arguments.
if test -z "$OUTPUT_SUFFIX" ; then
    OUTPUT_SUFFIX="mat"
fi

VERBOSE_OCTAVE=
if test "$VERBOSE_COUNT" -gt 2 ; then
    VERBOSE_OCTAVE="--verbose"
fi

# Now go through the list of parameters again to process the files. Assume that
# everything that does not match a valid option string is an input filename.
if test "$VERBOSE_COUNT" -gt 0 ; then
    echo "Converting:"
fi
for ARG ; do
    case $ARG in
        --help|-h|--output=*|-o=*|--suffix=*|-s=*|--force|-f|--verbose|-v)
            # Ignore option string
            ;;
        *)
            # Assume this argument is an input filename, so process:
            # Need to first create an output filename for the input file.
            INFILE="$ARG"
            if test -n "$OUTPUT_FILENAME" ; then
                OUTFILE="$OUTPUT_FILENAME"
            else
                OUTFILE="`dirname $INFILE`/`basename $INFILE | sed 's|\..*||'`.$OUTPUT_SUFFIX"
            fi
            if test "$VERBOSE_COUNT" -gt 0 ; then
                echo "$INFILE => $OUTFILE"
            fi
            if test -f "$OUTFILE" ; then
                if test "$FORCE_OVERWRITE" = no ; then
                    echo "WARNING: File '$OUTFILE' already exists so skipping conversion for file '$INFILE'." 1>&2
                    continue
                fi
            fi
            if test "$VERBOSE_COUNT" -gt 1 ; then
                octave $VERBOSE_OCTAVE --path "$MFILE_PATH" --eval "data = readDVHfile('$INFILE'); save '$OUTFILE' data;"
            else
                octave $VERBOSE_OCTAVE --path "$MFILE_PATH" --eval "data = readDVHfile('$INFILE'); save '$OUTFILE' data;" > /dev/null
            fi
            ;;
    esac
done
