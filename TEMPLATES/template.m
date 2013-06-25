%data = readDVHfile(filename)
%
%Reads in a standard DVH output text file given by the filename and returns the
%contents in the header and structures cell array.
%
%data - this contains a structure of the parsed DVH file. It will contain two
%       fields, 'header' and 'structures', structured as follows:
%  header - this contains a Matlab struct of field values that form the global
%           header of the DVH file.
%  structures - this is a cell array of Matlab structs that each contain fields
%               relevant to a biological structure (organ).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    Particle Therapy Project Bergen (PTPB) - tools and models for research in
%    cancer therapy using particle beams.
%
%    Copyright (C) 2013 Particle Therapy Group Bergen
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Author(s): name1 <blabla1@bla.no>, name2 <blabla2@bla.no>,... 
