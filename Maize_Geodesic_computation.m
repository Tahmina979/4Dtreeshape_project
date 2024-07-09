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

% We take samples from the saved spatiotemporally registered 4D plants
% If you want to start with the raw shapes you have to run the codes for
% spatial and temporal registration first and then give the output of those
% into here as input

ctree_1 = load('seq5s_alligned.mat');

num_within_seq1=numel(ctree_1.aligned_q_tree_5);

for i=1:num_within_seq1
    q_tree_1{i}=ctree_1.aligned_q_tree_5{1,i};
end

ctree_2 = load('seq2s_alligned.mat');
num_within_seq2=numel(ctree_2.aligned_q_tree_2);

for i=1:num_within_seq2
    q_tree_2{i}=ctree_2.aligned_q_tree_2{1,i};
end
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
   show_all_shapes_geod_maize(PC,count);
   
end

