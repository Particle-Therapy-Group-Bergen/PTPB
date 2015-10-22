#!/bin/sh

###############################################################################
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
###############################################################################

# Launches the make process with to run in parallel over as many CPUs as are
# available on this machine.

OCTAVE="octave -q --path ~/bin/PTPB_mfiles"

CORE_COUNT=`python -c "import multiprocessing; print multiprocessing.cpu_count();"`

if [ ! \( -d data -o -L data \) ] ; then
    echo "ERROR: the data directory has not been setup." 1>&2
    exit 1
fi

make -j "$CORE_COUNT" OCTAVE="$OCTAVE" $@
