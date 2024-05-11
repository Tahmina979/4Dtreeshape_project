clc;
clear;
close all;

addpath('utils_data','utils_draw','utils_funcs',"utils_statModels", "OpenCurvesRn");

lam_m = 1; 
lam_s = 1;
lam_p = 1;

gcf = figure;

set(gcf,'visible', 'on');
set(gcf, 'color', 'w');
view(0, 0);
axis equal; hold on;
box on; hold on;

set(gca,'XColor', 'none','YColor','none','ZColor','none');

data_path1= 'utils_data/NeuroData/maize01_seq/';
[all_qCompTrees1, all_compTrees1] = maize_load_from_mat_file(data_path1);
num_within_seq1=numel(all_qCompTrees1);
maize_show_all_shapes_sp_reg(all_compTrees1,1);


data_path2= 'utils_data/NeuroData/maize03_seq/'; 
[all_qCompTrees2, all_compTrees2] = maize_load_from_mat_file(data_path2);
num_within_seq2=numel(all_qCompTrees2);
maize_show_all_shapes_sp_reg(all_compTrees2,2)

fprintf('Please wait around 15-20 seconds to perform spatial registration within and accross growing plnats.....\n');

%registration within_sequence
tic
aligned_q_tree_1{1}=all_qCompTrees1{1};
for i=2:num_within_seq1
    Q1=all_qCompTrees1{i};
    Q2 = aligned_q_tree_1{i-1};
    [G,Q1p, Q2p] = ReparamPerm_qCompTrees_rad_4layers_v2(Q1, Q2, lam_m, lam_s, lam_p);
    aligned_q_tree_1{i}=Q1p;
end


j=1;
for i=1:num_within_seq2
    bol=isempty(all_qCompTrees2{i});
    if bol==1
        continue;
    end
   aligned_q_tree_2{j}=all_qCompTrees2{i};
   j=j+1;
end

if j==1
    num_within_seq2=numel(all_qCompTrees2);
    aligned_q_tree_2=all_qCompTrees2;
else
    num_within_seq2=numel(aligned_q_tree_2);
end

for i=2:num_within_seq2
    Q1=aligned_q_tree_2{i};
    Q2 = aligned_q_tree_2{i-1};
    [G,Q1p, Q2p] = ReparamPerm_qCompTrees_rad_4layers_v2(Q1, Q2, lam_m, lam_s, lam_p);
    aligned_q_tree_2{i}=Q1p;
end

%Equalize length
[Mu,eigenVector,EVals,tree_pca,used_qCompTrees_Ready]=transfer_to_PCA_geod_space(all_qCompTrees1,aligned_q_tree_2,num_within_seq1,num_within_seq2);

X1=tree_pca(1:num_within_seq1,:)';
X2=tree_pca(num_within_seq1+1:num_within_seq1+num_within_seq2,:)';

X2=ReSampleCurve(X2,num_within_seq1);

for i=1:size(X2,2)
       tree1=X2(:,i)'*eigenVector';

       tree=tree1+Mu';

       recon{i}=tree;

 end

for j=1:size(X2,2)

   aligned_q_tree_2{j}=unflattenCompTree_4layers_rad(recon{j}, lam_m, lam_s, lam_p, used_qCompTrees_Ready{1});
end

num_within_seq2=numel(aligned_q_tree_2);


%Registration across sequence
for i=1:num_within_seq1
    Q1=aligned_q_tree_2{i};
    Q2=aligned_q_tree_1{i};
    [G,Q1p, Q2p] = ReparamPerm_qCompTrees_rad_4layers_v2(Q1, Q2, lam_m, lam_s, lam_p);
    aligned_q_tree_2{i}=Q1p;
end

toc

for i=1:num_within_seq1
    PC1{i}=qCompTree_to_CompTree_rad_4layers(aligned_q_tree_1{i});
end

maize_show_all_shapes_sp_reg(PC1,3);

for i=1:num_within_seq2
    PC2{i}=qCompTree_to_CompTree_rad_4layers(aligned_q_tree_2{i});
end

maize_show_all_shapes_sp_reg(PC2,4);
