function [results, uncertainty] = makeResultsTable(organtech, M, N)
%[results, uncertainty] = makeResultsTable(organtech, M, N)
%
%Produces a table of results for a certain technique and organ combination.
%These are printed to standard output.
%
%organtech - This should be one of: 'cion_bladder', 'cion_rectum',
%            'proton_bladder', 'proton_rectum'
%M - The number of boot-strap samples to create to estimate the uncertainty of
%    the statistical parameters.
%N - The number of input samples to compute the parameters from.
%
%The following outputs are produced:
%
%results - The values of the statistical parameters, where each row contains the
%          results of calcDistribStats, with one row per loaded file.
%uncertainty - The uncertainties of the statistical parameters with the same
%              structure as for 'results'.

if ~ exist('organtech')
    organtech = 'cion_bladder';
end
if ~ exist('M')
    M = 100;
end

filepatterns = {
    '%s_sample_data.mat',
    '%s_sample_rbe_one_data.mat',
    '%s_sample_no_frac_data.mat',
    '%s_sample_no_alpha_data.mat',
    '%s_sample_no_beta_data.mat',
    '%s_sample_no_RBEmin_data.mat',
    '%s_sample_no_RBEmax_data.mat',
    '%s_sample_no_bin_spread_data.mat',
    '%s_sample_no_spread_data.mat'
};

for n = 1:length(filepatterns)
    filepattern = filepatterns{n};
    filename = sprintf(filepattern, organtech);
    if ~ exist('N')
        [R, U, D] = calcDistribStats(filename, M);
    else
        [R, U, D] = calcDistribStats(filename, M, N);
    end
    results(n,:) = R;
    uncertainty(n,:) = U;
    if n == 1
        header = {
                '2.5% quant',
                '    median',
                '97.5% quant',
                '      mean',
                ' std. dev.',
                'int. size.'
            };
        for k = 1:length(header)
            if k == 1
                printf('%s +/- std. dev.', header{k});
            else
                printf('\t%s +/- std. dev.', header{k});
            end
        end
        printf('\tfilename\n');
    end
    for k = 1:length(R)
        intsize = R(k);
        interr = U(k);
        if k == 1
            printf('%.9g +/- %.9g', intsize, interr);
        else
            printf('\t%.9g +/- %.9g', intsize, interr);
        end
    end
    printf('\t%s\n', filename);
end
