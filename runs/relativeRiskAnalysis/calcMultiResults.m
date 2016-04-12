function [results, uncertainty, datacube] = calcMultiResults(M, N)

if ~ exist('M')
    M = 100;
end

filenames = {
        'cion_bladder_sample_multiA_RBE_data.mat',
        'cion_rectum_sample_multiA_RBE_data.mat',
        'proton_bladder_sample_multiA_RBE_data.mat',
        'proton_rectum_sample_multiA_RBE_data.mat',
        'cion_bladder_sample_multiB_RBE_data.mat',
        'cion_rectum_sample_multiB_RBE_data.mat',
        'proton_bladder_sample_multiB_RBE_data.mat',
        'proton_rectum_sample_multiB_RBE_data.mat'
    };

for n = 1:length(filenames)
    filename = filenames{n};
    if ~ exist('N')
        [R, U, D] = calcDistribStats(filename, M);
    else
        [R, U, D] = calcDistribStats(filename, M, N);
    end
    results(n,:) = R;
    uncertainty(n,:) = U;
    datacube(n,:,:) = D;
    printf('%g +/- %g .. %g +/- %g .. %g +/- %g  %s\n',
           results(n,1), uncertainty(n,1),
           results(n,4), uncertainty(n,4),
           results(n,3), uncertainty(n,3),
           filename);
end
