function [alpha, beta, RBEmin, RBEmax, RR, RRerr] = makeRelativeRiskHyperCube(
                    organ, patients, filepat1, filepat2, n1, n2, scale1, scale2,
                    RBEmin_vec, RBEmax_vec)


% Input DVH data file pattern 1.
if ~ exist('filepat1')
    filepat1 = 'data/VMATdvh/vmat%d.mat';
end

% number of fractions 1.
if ~ exist('n1')
    n1 = 25;
end

% Dose scaling factor for DVH 1. Dose will be multiplied by this.
if ~ exist('scale1')
    scale1 = 1;
end

% Input DVH data file pattern 2.
if ~ exist('filepat2')
    filepat2 = 'data/CionDataPhysicalDose/HUH%dphysical_dvh.mat';
end

% number of fractions 2.
if ~ exist('n2')
    n2 = 12;
end

% Dose scaling factor for DVH 2. Dose will be multiplied by this.
if ~ exist('scale2')
    scale2 = 1;
end

if ~ exist('patients')
    patients = [12, 33, 35, 36, 37, 39, 41, 42, 43, 44];
end

if ~ exist('organ')
    organ = 'Bladder';
end

if ~ exist('RBEmin_vec')
    RBEmin_vec = [1.05 1.15 1.25 1.35];
end

if ~ exist('RBEmax_vec')
    RBEmax_vec = [2 4 6 8];
end

opts = struct('integration_method', 'trapz',
              'integration_tolerance', 1e-4,
              'interpolation_method', 'pchip');
namemap = {
        'Bladder_P', 'Bladder';
        'Rectum_P_MT', 'Rectum';
    };

% Load DVH data from files.
dvh1list = {};
dvh2list = {};
for k = 1:length(patients)
    % Relative Risk = file1/file2
    file1 = sprintf(filepat1, patients(k));
    file2 = sprintf(filepat2, patients(k));
    dvh1 = getDoseVolumeHistogram(file1, namemap, organ).(organ).datapoints;
    dvh1(:,1) = dvh1(:,1) .* scale1;
    dvh1list{k} = dvh1;
    dvh2 = getDoseVolumeHistogram(file2, namemap, organ).(organ).datapoints;
    dvh2(:,1) = dvh2(:,1) .* scale2;
    dvh2list{k} = dvh2;
end

alpha_vec = 0.01:0.02:0.8;
beta_vec = 0.010:0.001:0.050;

alpha = zeros(length(alpha_vec), length(beta_vec),
              length(RBEmin_vec), length(RBEmax_vec));
beta = alpha;
RBEmin = alpha;
RBEmax = alpha;
RR = alpha;
RRerr = alpha;

for x = 1:length(alpha_vec)
for y = 1:length(beta_vec)
for z = 1:length(RBEmin_vec)
for w = 1:length(RBEmax_vec)
    alpha(x, y, z, w) = alpha_vec(x);
    beta(x, y, z, w) = beta_vec(y);
    RBEmin(x, y, z, w) = RBEmin_vec(z);
    RBEmax(x, y, z, w) = RBEmax_vec(w);
    [rr, rrerr] = CalculateMeanRR(dvh1list, dvh2list, opts, n1, n2,
                                  alpha(x, y, z, w), beta(x, y, z, w),
                                  RBEmin(x, y, z, w), RBEmax(x, y, z, w));
    RR(x, y, z, w) = rr;
    RRerr(x, y, z, w) = rrerr;
end
end
end
end

return;


function [rr, rrerr] = CalculateMeanRR(dvh1list, dvh2list, opts, n1, n2,
                                       alpha, beta, RBEmin, RBEmax)
RR = zeros(length(dvh1list), 1);
for k = 1:length(RR)
    % We have to make sure we tie the distribution end points to 1 and 0 for
    % the volume ratios, otherwise the integration becomes unstable.
    dvh1 = dvh1list{k};
    dvh1 = [0, 1; dvh1; dvh1(length(dvh1),1)+1, 0];
    dvh2 = dvh2list{k};
    dvh2 = [0, 1; dvh2; dvh2(length(dvh2),1)+1, 0];
    RR(k) = RelativeRisk('LinearQuad', dvh1, dvh2, opts, n1, n2, alpha, beta,
                         RBEmin, RBEmax);
end
rr = mean(RR);
rrerr = std(RR)/sqrt(length(RR));
fprintf(stdout,
        'relative risk = %g +/- %g\talpha = %g\tbeta = %g\tRBEmin = %g\tRBEmax = %g\n',
        rr, rrerr, alpha, beta, RBEmin, RBEmax);
fflush(stdout);
return;
