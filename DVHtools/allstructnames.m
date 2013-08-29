% Writes all structure names in DVH file to screen

% From within MatLab/Octave console run:
% allstructnames('InputDVHfile')


function allstructnames(filename) 

% printing all organ names

% reads DVHinput and prints patient ID
fprintf('Reading file "%s"...\n', filename);
fflush(stdout);
data = readDVHfile(filename);
fprintf('Done.\n');
fflush(stdout);
fprintf('DVH for patient: %s\n', data.header.patientName);

for j=1:length(data.structures)
	structure=data.structures{j}.structName;
	fprintf('number %d \n', j);
	fprintf('structure %s \n', structure);
end


