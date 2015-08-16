#!/usr/bin/env python

# Utility program to quickly convert XiO DVH text output to Varian output format.

from __future__ import print_function
import sys
import os
import argparse


def Convert(infile, outfile):
    print("Converting " + infile.name + " => " + outfile.name)
    lines = infile.read().split('\n')
    state = 'start'
    lineno = 0
    dvh_count = 0
    data = []
    patient_id = None
    patient_name = None
    organ = None
    resolution = None
    bins = None
    date = None
    datarows = []
    # Parse input:
    for line in lines:
        lineno += 1
        if not line:  # skip empty lines.
            continue
        elif "XiO Dose Volume Histogram" == line:
            if state != 'start':
                elem = {'organ': organ, 'resolution': resolution, 'bins': bins,
                        'datarows': datarows}
                data.append(elem)
            state = 'new_dvh'
            dvh_count += 1
            organ = None
            resolution = None
            bins = None
            datarows = []
        elif state == 'new_dvh':
            state = 'got_id'
            if patient_id is None:
                patient_id = line
            else:
                if line != patient_id:
                    raise Exception("Error: Patient ID different in " + str(infile.name) \
                                    + ":" + str(lineno) + ", found: " + line)
        elif state == 'got_id':
            state = 'got_name'
            patient_name = line
            if patient_name is None:
                patient_name = line
            else:
                if line != patient_name:
                    raise Exception("Error: Patient name different in " + str(infile.name) \
                                    + ":" + str(lineno) + ", found: " + line)
        elif state == 'got_name':
            state = 'got_organ'
            organ = line
        elif state == 'got_organ':
            state = 'got_res'
            fields = line.split('\t')
            if len(fields) == 2 and fields[1] == "mm resolution":
                resolution = fields[0]
            else:
                raise Exception("Error: Resolution field format is not correct in " \
                                + str(infile.name) + ":" + str(lineno) + ", found: " + line)
        elif state == 'got_res':
            state = 'got_bins'
            fields = line.split('\t')
            if len(fields) == 2 and fields[1] == "bins":
                bins = fields[0]
            else:
                raise Exception("Error: Bins field format is not correct in " \
                                + str(infile.name) + ":" + str(lineno) + ", found: " + line)
        elif state == 'got_bins':
            state = 'got_date'
            date = line
            if date is None:
                date = line
            else:
                if line != date:
                    raise Exception("Error: Date different in " + str(infile.name) \
                                    + ":" + str(lineno) + ", found: " + line)
        elif state == 'got_date' and "Min. Bin Dose (cGy),  Bin Volume (cc)" in line:
            state = 'data_section'
        elif state == 'data_section':
            fields = line.split(',')
            if len(fields) == 2:
                s = fields[0].strip()
                dose = float(s)
                dose /= 100
                s = fields[1].strip()
                volume = float(s)
                datarows.append( [dose, volume] )
            else:
                print("WARNING: unexpected field format " + str(infile.name) \
                      + ":" + str(lineno) + ": " + str(line), file=sys.stderr)
        else:
            print("WARNING: unexpected line in " + str(infile.name) + ":" \
                  + str(lineno) + ": " + str(line), file=sys.stderr)
    # Make sure we add the last organ in the file when we reach the end of the file.
    if state == 'data_section' and organ is not None and resolution is not None \
       and bins is not None:
        elem = {'organ': organ, 'resolution': resolution, 'bins': bins,
                'datarows': datarows}
        data.append(elem)
    # Generate output:
    if patient_name is None or patient_id is None or date is None:
        raise Exception("Error: Could not find patient name, ID or date fields in "\
                        + infile.name + ".")
    print("Patient Name :", patient_name, file=outfile)
    print("Patient ID   :", patient_id, file=outfile)
    print("Comment      :", "coverted from XiO", file=outfile)
    print("Date         :", date, file=outfile)
    print("Exported by  :", "DVH_tools", file=outfile)
    print("Type         :", "Cumulative Dose Volume Histogram",
            file=outfile)
    print("Description  :", "This is a converted file.",
            file=outfile)
    print("", file=outfile)
    print("Plan sum     :", "unknown", file=outfile)
    print("Prescribed dose [Gy]:", "not defined", file=outfile)
    print("% for dose (%):", "not defined", file=outfile)
    for elem in data:
        total_volume = 0
        for row in elem['datarows']:
            total_volume += row[1]
        print("", file=outfile)
        print("Structure               :", elem['organ'], file=outfile)
        print("Approval Status         :", "Unapproved", file=outfile)
        print("Plan                    :", "unknown", file=outfile)
        print("Course                  :", "unknown", file=outfile)
        print("Volume [cm^3]           :", total_volume, file=outfile)
        print("Dose Cover.[%]          :", "100.0", file=outfile)
        print("Sampling Cover.[%]      :", "100.0", file=outfile)
        print("Min Dose [Gy]           :", "N/A", file=outfile)
        print("Max Dose [Gy]           :", "N/A", file=outfile)
        print("Mean Dose [Gy]          :", "N/A", file=outfile)
        print("Modal Dose [Gy]         :", "N/A", file=outfile)
        print("Median Dose [Gy]        :", "N/A", file=outfile)
        print("STD [Gy]                :", "N/A", file=outfile)
        print("Equiv. Sphere Diam. [cm]:", "N/A", file=outfile)
        print("Conformity Index        :", "N/A", file=outfile)
        print("Gradient Measure [cm]   :", "N/A", file=outfile)
        #print("Resolution [mm]         :", elem['resolution'], file=outfile)
        #print("Bins                    :", elem['bins'], file=outfile)
        print("", file=outfile)
        header = "{0:>16} {1:>40}".format("Dose [Gy]", "Ratio of Total Structure Volume [%]")
        print(header, file=outfile)
        volume = total_volume
        for row in elem['datarows']:
            volume_ratio = volume / total_volume * 100
            line = "{0:>16f} {1:>40f}".format(row[0], volume_ratio)
            print(line, file=outfile)
            volume -= row[1]  # increment at end to follow how XiO does it.


# Setup the command line argument parser with help messages:
argparser = argparse.ArgumentParser(description="Converts XiO DVH text file"
    " output to a format used by Varian.")
argparser.add_argument("infile", metavar="<input file>", nargs='+', type=open,
    help="One or more input files to convert.")
argparser.add_argument( "-o", "--output", dest="outfilename",
    metavar="<filename>", default=None,
    help="The name of the output file. Can only be used with a single input"
    " file specified.")
argparser.add_argument( "-s", "--suffix", dest="file_suffix",
    metavar="<ext>", default=None,
    help="File name suffix to use for the output file.")

# Parse the command line arguments:
try:
    args = argparser.parse_args()
    if len(args.infile) > 1 and not (args.outfilename is None):
        raise Exception("Error: If the -o | --output option is used then only"
                        " one input file can be given.")
    for infile in args.infile:
        if args.outfilename is None:
            basename, ext = os.path.splitext(infile.name)
            if args.file_suffix is None:
                new_ext = "xio"
            else:
                new_ext = args.file_suffix
            outfilename = basename + '.' + new_ext
        else:
            outfilename = args.outfilename
        outfile = open(outfilename, "w")
        Convert(infile, outfile)
        outfile.close()
except Exception as e:
    print(str(e), file=sys.stderr)
