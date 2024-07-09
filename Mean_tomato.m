clc;
clear;
%close all;

addpath('utils_data','utils_draw',"OpenCurvesRn",'utils_funcs',"utils_statModels",'Registered_tomato_plants');

lam_m = 1; 
lam_s = 1;
lam_p = 1;

gcf = figure;%('position', [0, 0, 1500, 800]);

set(gcf,'visible', 'on');
set(gcf, 'color', 'w');
view(0, 0);
axis equal; hold on;
box on; hold on;

set(gca,'XColor', 'none','YColor','none','ZColor','none');

% We take all the spatiotemporally registered 4D shapes. Use the codes in
% "registration_code" folder to perform the spatiotemporal registration on
% a set of 4D tree shapes 

ctree_1 = load('Registered_tomato_plants/seq1s_alligned.mat');

ctree_2 = load('Registered_tomato_plants/seq2s_alligned.mat');

ctree_3 = load('Registered_tomato_plants/seq3s_alligned.mat');

ctree_4 = load('Registered_tomato_plants/seq4s_alligned.mat');

ctree_5 = load('Registered_tomato_plants/seq5s_alligned.mat');

ctree_6 = load('Registered_tomato_plants/seq6s_alligned.mat');

ctree_7 = load('Registered_tomato_plants/seq7s_alligned.mat');

num_within_seq1=numel(ctree_1.aligned_q_tree_1);
num_within_seq2=numel(ctree_2.aligned_q_tree_2);
num_within_seq3=numel(ctree_3.aligned_q_tree_3);
num_within_seq4=numel(ctree_4.aligned_q_tree_4);
num_within_seq5=numel(ctree_5.aligned_q_tree_5);
num_within_seq6=numel(ctree_6.aligned_q_tree_6);
num_within_seq7=numel(ctree_7.aligned_q_tree_7);



[Mu,eigenVector,EVals,tree_pca,used_qCompTrees_Ready]=transfer_to_PCA_space(ctree_1,ctree_2,ctree_3,ctree_4,ctree_5,ctree_6,ctree_7,num_within_seq1,num_within_seq2,num_within_seq3,num_within_seq4,num_within_seq5,num_within_seq6,num_within_seq7);


%reshaping curves
start=1;
for i=1:7 % 7 is the total number of 4D sequence in the set
    ending=start+9;
    X(:,i)=reshape(tree_pca(start:ending,:)',[],1);
    start=ending+1;
end


%PCA on curves

[Mu_s, eigenVectors_s, EVals_s] = performEigenAnalysis(X);

% Mu_s is the mean curve

Mean=reshape(Mu_s,[size(tree_pca,2),num_within_seq1]) %all the numbers are same, we can use any
for i=1:num_within_seq1
       tree1=Mean(:,i)'*eigenVector'; 
       tree=tree1+Mu';
       recon{i}=tree;



PC1{i}=unflattenCompTree_4layers_rad(recon{i}, lam_m, lam_s, lam_p, used_qCompTrees_Ready{1});
PC{i}=qCompTree_to_CompTree_rad_4layers(PC1{i});

end
show_all_shapes_modes_tomato(PC,1);



return;
