clc;
clear;
close all


%% ############ DATA ############

% Botanical Trees
addpath('utils_data');
addpath('utils_data/utils_data_botanTrees');
addpath('utils_statModels')

data_path = 'botanTrees_txtskl_SGP18';              % botanTrees_txtskl_SGP18
[all_qCompTrees, all_compTrees] = load_botanTrees_rad(data_path);

% Hard coded
G_1_idxes = linspace(1, 6, 6);
G_2_idxes = linspace(7, 9, 3);
G_3_idxes = linspace(10, 32, 23);
G_4_idxes = linspace(33, 38, 6);

% used_idxes = G_4_idxes;
meanTree_G1_idxes = [1, 7, 10];
meamTree_G2_idxes = [2, 6, 29];
meanTree_G3_idxes = [15, 17, 24, 30, 31];

meanTree_G4_idxes = [3, 4, 5];
meanTree_G5_idxes = [33, 34, 35, 36, 37];

% [used_qCompTrees, used_compTrees] = augment_botanTrees_rad(all_qCompTrees(used_idxes), ...
%                                                            all_compTrees(used_idxes) );

used_idxes = meanTree_G5_idxes;

used_qCompTrees = all_qCompTrees(used_idxes);
used_compTrees = all_compTrees(used_idxes);
tNum = length(used_qCompTrees);
% 
addpath('utils_draw')
run showAll_compBotanTrees_ownMade4layers.m
% return

% % Load augmented trees
% load('rand_sampleQ')
% used_qCompTrees = [used_qCompTrees, rand_sampleQ];
% tNum = length(used_qCompTrees);

%% ######### Statistical Models #########

% %%%%%% Mean Computation %%%%%%
addpath('utils_statModels')

% parameters
lam_m = 0.1; 
lam_s = 0.1;
lam_p = 1;
Nitr = 3;

fprintf('Q1 to Q2 (with perm), lam_m:%.2f, lam_s:%.2f, lam_p:%.2f\n', ...
                                                    lam_m, lam_s, lam_p);


% === Loops: Pad, Align tree pairs; Compuate geodesic ====
tm_all = tic;
qMean = used_qCompTrees{1};
for i = 2: tNum
    
    tm1 = tic;
    Q1 = qMean;
    Q2 = used_qCompTrees{i};
    
    % ---- Pad and Align trees, NO procrustes analysis ---
    [G,Q1p, Q2p] = ReparamPerm_qCompTrees_rad_4layers_v2(Q1, Q2, lam_m, lam_s, lam_p);
    
    T1 = toc(tm1); fprintf('Loop %d: Pad and Align trees - done, timecost:%.4f secs\n', (i-1), T1);
    
    % --- Align two trees and compute correspondence ---
    % [Q1p, G, ~, Q2p] = AlignFull_qComplexTree_4layers(Q1, Q2, lam_m,lam_s,lam_p, Nitr); 

    % --- Compute Geodesic ---
    tm2 = tic;
    stp1 = i+1;
    [A10, qA10] = GeodComplexTreesPrespace_rad_4layers(Q1p, Q2p, stp1);
    qMean = qA10{2};
    
    T2 = toc(tm2); fprintf('Loop %d: Geodesic computation - done, timecost:%.4f secs\n',(i-1), T2);

end

T3 = toc(tm_all); fprintf('Mean Loops - done, timecost:%.4f secs\n', T3);

% --- Visualization, A10{2} is final mean ----------
addpath('utils_draw')
run showInputAndMean_compBotanTrees_4layers.m
return;

% --- Save Objs --------
addpath('GetOBJ')
Data = [used_compTrees(1:end), A10{2}];
obj_folder = saveInputAndMeanObjs_compTrees_rad_4layers(Data, [data_path,'--'], used_idxes);


% ===== Render Objs =====
addpath('RenderOBJ')
addpath('RenderOBJ/func_render/');
addpath('RenderOBJ/func_other/');
run renderBotanInputAndMeanObjs.m



return;


