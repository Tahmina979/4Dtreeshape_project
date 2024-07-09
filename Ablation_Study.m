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

%target
data_path2= 'utils_data/NeuroData/tom_3_seq_altered/';
[all_qCompTrees2, all_compTrees2] = load_from_mat_file(data_path2);

num_within_seq2=numel(all_compTrees2);

show_all_shapes_reg(all_compTrees2,2);
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
aligned_q_tree_2{1}=aligned_q_tree_2{1};
%registration within_sequence
for i=2:num_within_seq2
   Q1=aligned_q_tree_2{i};
  Q2 = aligned_q_tree_2{i-1};
    [G,Q1p, Q2p] = ReparamPerm_qCompTrees_rad_4layers_v2(Q1, Q2, lam_m, lam_s, lam_p);
    aligned_q_tree_2{i}=Q1p;
   
end

%Taking sequences registered spatially within them,SourceGT
ctree_1 = load('seq3_alligned.mat');

num_within_seq1=numel(ctree_1.aligned_q_tree_3);

for i=1:num_within_seq1
    aligned_q_tree_1{i}=ctree_1.aligned_q_tree_3{1,i};
   display_source{i}=qCompTree_to_CompTree_rad_4layers(aligned_q_tree_1{i});
end
show_all_shapes_reg(display_source,1);


%test in q_space
[Muq,eigenVectorq,EValsq,tree_pcaq,used_CompTrees_Readyq,qXq]=transfer_to_PCA_geod_space(aligned_q_tree_1,aligned_q_tree_2,num_within_seq1,num_within_seq2);

X1q=tree_pcaq(1:num_within_seq1,:)';
X2q=tree_pcaq(num_within_seq1+1:num_within_seq1+num_within_seq2,:)';
X2q=ReSampleCurve(X2q,num_within_seq1);

for i=1:size(X2q,2)
       tree1=X2q(:,i)'*eigenVectorq';

       tree=tree1+Muq';

       recon{i}=tree;

 end

for j=1:size(X2q,2)

  aligned_q_tree_2{j}=unflattenCompTree_4layers_rad(recon{j}, lam_m, lam_s, lam_p, used_CompTrees_Readyq{1});
end
num_within_seq2=numel(aligned_q_tree_2);

%registration_accross_seq
for i=1:num_within_seq1
    
    Q1=aligned_q_tree_2{i};
  Q2=aligned_q_tree_1{i};
   
  [G,Q1p, Q2p] = ReparamPerm_qCompTrees_rad_4layers_v2(Q1, Q2, lam_m, lam_s, lam_p);
  aligned_q_tree_2{i}=Q1p;
end

[Muq,eigenVectorq,EValsq,tree_pcaq,used_CompTrees_Readyq,qXq]=transfer_to_PCA_geod_space(aligned_q_tree_1,aligned_q_tree_2,num_within_seq1,num_within_seq2)
X1q=tree_pcaq(1:num_within_seq1,:)';
X2q=tree_pcaq(num_within_seq1+1:num_within_seq1+num_within_seq2,:)';
[distq,X2nq,q2n,X1q,X2q,q1, q2,gamIq, distbeforeq]=mygeod(X1q, X2q);

for i=1:size(X2nq,2)
     tree1=X2nq(:,i)'*eigenVectorq';

       tree=tree1+Muq';

       recon{i}=tree;

 end


for j=1:size(X2nq,2)
  in_curve_q{j}=unflattenCompTree_4layers_rad(recon{j}, lam_m, lam_s, lam_p, used_CompTrees_Readyq{1});
  in_curve_q_space{j}=qCompTree_to_CompTree_rad_4layers(in_curve_q{j});
end
show_all_shapes_reg(in_curve_q_space,3);

%test in original space

%warp_target
i=1;
for j=1:numel(aligned_q_tree_2)


  if j==2||j==4||j==9
      continue;
  end
  aligned_q_tree_2_altered{i}=aligned_q_tree_2{j};
  i=i+1;
 
end

num_within_seq2=numel(aligned_q_tree_2_altered);

[Mu,eigenVector,EVals,tree_pca,used_CompTrees_Ready,qX]=transfer_to_PCA_geod_space_ablation(aligned_q_tree_1,aligned_q_tree_2_altered,num_within_seq1,num_within_seq2);
X1=tree_pca(1:num_within_seq1,:)';
X2=tree_pca(num_within_seq1+1:num_within_seq1+num_within_seq2,:)';
X2=ReSampleCurve(X2,num_within_seq1);
num_within_seq2=num_within_seq1;
[dist,X2n,X1,X2,gamI, distbefore]=mygeod_org(X1, X2);
for i=1:size(X2n,2)
     tree1=X2n(:,i)'*eigenVector';

       tree=tree1+Mu';

       recon{i}=tree;

end
for j=1:size(X2n,2)

  in_org{j}=unflattenphenoTree_4layers_rad(recon{j}, lam_m, lam_s, lam_p, used_CompTrees_Ready{1});

end
show_all_shapes_reg(in_org,4);





