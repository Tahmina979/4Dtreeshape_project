
fp1 = fopen('prog_info.txt', 'w');

% load qX
start_t1 = tic;
load('qX.mat');
T1= toc(start_t1);
fprintf(fp1, 'qX loaded-Done..., timecost: %f s\n', T1);

% covariance matrix
start_t2 = tic;
covMat = cov(qX);

T2 = toc(start_t2);
fprintf(fp1, 'covMat computation-Done..., timecost: %f s\n', T2);

qX_mean = mean(qX);

% SVD
start_t3 = tic;
% [V,E,U] = svd(covMat);
[V,E,U] = svds(covMat, 10);
save('V_E_U.mat', 'V', 'E', 'U');

T3 = toc(start_t3);
fprintf(fp1, 'SVD computation-Done..., timecost: %f s', T3);