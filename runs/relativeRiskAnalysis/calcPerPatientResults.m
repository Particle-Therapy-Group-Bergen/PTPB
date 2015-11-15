function [results, uncertainty, datacube] = calcPerPatientResults(patientIDs, filepattern, M, N)
%[results, uncertainty, datacube] = calcPerPatientResults(patientIDs, filepattern, M, N)
%
%Calculates various statistical parameters for the sample data for each patient.
%
%patientIDs - A vector of patient IDs to process.
%filepattern - A printf compatible file name pattern to load sample data from.
%              The pattern will be filled with each entry from 'patientIDs'.
%M - The number of boot-strap samples to create to estimate the uncertainty of
%    the statistical parameters.
%N - The number of input samples to compute the parameters from.
%
%The following outputs are produced:
%
%results - The values of the statistical parameters, where each row contains the
%          results of calcDistribStats, with one row per patient.
%uncertainty - The uncertainties of the statistical parameters with the same
%              structure as for 'results'.
%datacube - The statistical parameters samples used to estimate the
%           uncertainties of the final results. There is one matrix per patient,
%           which is selected in the first dimention.

if ~ exist('patientIDs')
    patientIDs = [12 33 35 36 37 39 41 42 43 44];
end
if ~ exist('filepattern')
    filepattern = 'patient_%d_cion_bladder_sample_data.mat';
end
if ~ exist('M')
    M = 100;
end

for n = 1:length(patientIDs)
    filename = sprintf(filepattern, patientIDs(n));
    [R, U, D] = calcDistribStats(filename, M, N);
    results(n,:) = R;
    uncertainty(n,:) = U;
    datacube(n,:,:) = D;
    printf('Patient %d: %g +/- %g .. %g +/- %g .. %g +/- %g\n', patientIDs(n),
           results(n,1), uncertainty(n,1),
           results(n,4), uncertainty(n,4),
           results(n,3), uncertainty(n,3));
end
