clc;
clear;
close all


%% ############ DATA ############
% NeuroTrees
addpath('utils_data', 'utils_draw');

data_path = 'NeuroData/chen/CNG version/';
% data_path = 'NeuroData/wu/CNG version/';
[all_qCompTrees, all_compTrees] = load_neuroTrees_rad(data_path);

%% ######### Statistical Models #########

% ##### Geodesic Computation #####
addpath('utils_statModels')
idx1 =22; idx2 = 60;

Q1 = all_qCompTrees{idx1};
Q2 = all_qCompTrees{idx2};

% parameters
lam_m = 0.2; 
lam_s = 1;
lam_p = 0.2;

fprintf('Q1 to Q2 (with perm), lam_m:%.2f, lam_s:%.2f, lam_p:%.2f\n', ...
                                                    lam_m, lam_s, lam_p);


% ===== Pad and Align trees =====
tm1 = tic;
[G,Q1p, Q2p] = ReparamPerm_qCompTrees_rad_4layers_v2(Q1, Q2, lam_m, lam_s, lam_p);

T1 = toc(tm1); fprintf('Pad and Align trees - done, timecost:%.4f secs\n', T1); 


% ===== Compute the geodesic，map back to Euclidean space =====
stp1 = 9;
tm2 = tic;
[A10, qA10] = GeodComplexTreesPrespace_rad_4layers(Q1p,Q2p,stp1);

T2 = toc(tm2); fprintf('Geodesic computation - done, timecost:%.4f secs\n', T2);

% % ===== Save Objs =====
addpath('GetOBJ')
obj_folder = saveGeoObjs_compTrees_rad_4layers(A10, data_path, idx1, idx2);

% ===== Render Objs =====
addpath('RenderOBJ')
addpath('RenderOBJ/func_render/');
addpath('RenderOBJ/func_other/');
run renderNeuroGeodObjs.m






