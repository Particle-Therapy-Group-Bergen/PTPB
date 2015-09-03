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

# Top level make file to produce all calculations and plots.

OCTAVE = octave -q --path ~/bin/PTPB_mfiles

.PHONY: all clean cleanall

all:
	$(MAKE) -f hypercubes.makefile OCTAVE="$(OCTAVE)"
	$(MAKE) -f scan_plots.makefile OCTAVE="$(OCTAVE)"

clean:
	$(MAKE) -f hypercubes.makefile clean
	$(MAKE) -f scan_plots.makefile clean

cleanall:
	$(MAKE) -f hypercubes.makefile cleanall
	$(MAKE) -f scan_plots.makefile cleanall
