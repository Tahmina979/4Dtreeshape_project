clc;
clear;
close all;

addpath('utils_data','utils_draw',"OpenCurvesRn",'utils_funcs',"utils_statModels",'Registered_tomato_plants');


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

data_path1 = 'dataset/tom_5_seq/';
[all_qCompTrees1, all_compTrees1] = load_from_mat_file(data_path1);
num_within_seq1=numel(all_compTrees1);
show_shapes(all_compTrees1,1);


data_path2= 'dataset/tom_1_seq_altered/';
[all_qCompTrees2, all_compTrees2] = load_from_mat_file(data_path2);
num_within_seq2=numel(all_compTrees2);
show_shapes(all_compTrees2,2);

j=1;
for i=1:num_within_seq2
    bol=isempty(all_qCompTrees2{i});
    if bol==1
        continue;
    end
   aligned_q_tree_2{j}=all_qCompTrees2{i};
   j=j+1; 
end

num_within_seq2=numel(aligned_q_tree_2);

%registration within_sequence
for i=2:num_within_seq2
    Q1=aligned_q_tree_2{i};
    Q2 = aligned_q_tree_2{i-1};
    [G,Q1p, Q2p] = ReparamPerm_qCompTrees_rad_4layers_v2(Q1, Q2, lam_m, lam_s, lam_p);
    aligned_q_tree_2{i}=Q1p;
   
end
[Mu,eigenVector,EVals,tree_pca,used_qCompTrees_Ready,qX]=transfer_to_PCA_geod_space(all_qCompTrees1,aligned_q_tree_2,num_within_seq1,num_within_seq2);

X1=tree_pca(1:num_within_seq1,:)';
X2=tree_pca(num_within_seq1+1:num_within_seq1+num_within_seq2,:)';
tic
X2=ReSampleCurve(X2,num_within_seq1);

for i=1:size(X2,2)
       tree1=X2(:,i)'*eigenVector';

       tree=tree1+Mu';

       recon{i}=tree;

 end

for j=1:size(X2,2)

   aligned_q_tree_2{j}=unflattenCompTree_4layers_rad(recon{j}, lam_m, lam_s, lam_p, used_qCompTrees_Ready{1});
   PC{j}=qCompTree_to_CompTree_rad_4layers(aligned_q_tree_2{j});

end

show_shapes(PC,3);

num_within_seq2=size(X2,2);


[dist,X2n,q2n,X1,X2,q1, gamI, distbefore]=mygeod(X1, X2);
toc

for i=1:size(X2n,2)
      tree1=X2n(:,i)'*eigenVector';

      tree=tree1+Mu';

      recon{i}=tree;

 end

for j=1:size(X2n,2)

   aligned_q_tree_2{j}=unflattenCompTree_4layers_rad(recon{j}, lam_m, lam_s, lam_p, used_qCompTrees_Ready{1});
   PC2{j}=qCompTree_to_CompTree_rad_4layers(aligned_q_tree_2{j});

end

show_shapes(PC2,4);

data_path2= 'dataset/tom_1_seq/';
[all_qCompTrees2, all_compTrees2] = load_from_mat_file(data_path2);
num_within_seq2=numel(all_compTrees2);
show_shapes(all_compTrees2,5);