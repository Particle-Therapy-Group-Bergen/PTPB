###############################################################################
#
#    Particle Therapy Project Bergen (PTPB) - tools and models for research in
#    cancer therapy using particle beams.
#
#    Copyright (C) 2013 Particle Therapy Group Bergen
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

# This file is used to install the scripts and Matlab files into the users
# ~/bin directory by default. If PREFIX_DIR is given then this is used instead
# as the top level installation directory.

.PHONY: all install

ifeq ($(PREFIX_DIR),)
PREFIX_DIR=~/bin
endif

INSTALL_PATHS = DVHtools models scripts src

MATLAB_FILES=$(shell find $(INSTALL_PATHS) -type f -name "*.m")
SCRIPT_FILES=$(shell find $(INSTALL_PATHS) -type f -name "*.sh") \
             $(shell find $(INSTALL_PATHS) -type f -name "*.py")

all:

install:
	mkdir -p $(PREFIX_DIR)/PTPB_mfiles/
	install $(SCRIPT_FILES)  $(PREFIX_DIR)/
	for SCRIPT in $(notdir $(SCRIPT_FILES)) ; do \
		sed -i -e 's|@@MFILE_PATH@@|$(PREFIX_DIR)/PTPB_mfiles|' $(PREFIX_DIR)/$$SCRIPT ; \
	done
	install $(MATLAB_FILES) $(PREFIX_DIR)/PTPB_mfiles/
	@echo "INFO: Consider adding the following line to your ~/.bashrc for more convenience when using the Matlab scripts:"
	@echo "      alias octave='octave --path $(PREFIX_DIR)/PTPB_mfiles'"
