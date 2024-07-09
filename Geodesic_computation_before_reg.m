clc;
clear;
close all;

addpath('utils_data','utils_draw',"OpenCurvesRn",'utils_funcs',"utils_statModels");

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

% we compute geodesic on the raw data, if unequal length, we take interpolated one to make equal
% lenth 4D plant

ctree_1 = load('seq1_inter_24.mat');

num_within_seq1=numel(ctree_1.aligned_q_tree_1);

for i=1:num_within_seq1
    q_tree_1{i}=ctree_1.aligned_q_tree_1{1,i};
end

data_path1= 'utils_data/tom_3_seq/';
[q_tree_2, q_ctree_2] = load_from_mat_file(data_path1);
num_within_seq2=numel(q_tree_2);

[Mu,eigenVector,EVals,tree_pca,used_qCompTrees_Ready]=transfer_to_PCA_geod_space(q_tree_1,q_tree_2,num_within_seq1,num_within_seq2);

X1=tree_pca(1:num_within_seq1,:)';
X2=tree_pca(num_within_seq1+1:num_within_seq1+num_within_seq2,:)';


%geod computation
tic
Xgeod = computeGeodesic(X1,X2);
toc
count=0;
for k=1:size(Xgeod,1)
   for i=1:size(X2,2)
       tree1=Xgeod(k,:,i)*eigenVector';
       tree=tree1+Mu';
       recon{i}=tree;
   end

for j=1:size(recon,2)

   PC1{j}=unflattenCompTree_4layers_rad(recon{j}, lam_m, lam_s, lam_p, used_qCompTrees_Ready{1});
   PC{j}=qCompTree_to_CompTree_rad_4layers(PC1{j});

end
   count=count+1;
   show_all_shapes_geod(PC,count);
   
end

