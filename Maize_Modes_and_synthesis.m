clc;
clear;
%close all;

addpath('utils_data','utils_draw',"OpenCurvesRn",'utils_funcs',"utils_statModels",'Registered_maize_plants');

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




data_path1 = 'utils_data/NeuroData/maize01_seq_alt/';
[all_qCompTrees1, all_compTrees1] = maize_load_from_mat_file(data_path1);
num_within_seq1=numel(all_qCompTrees1);
maize_show_all_shapes_sp_reg(all_compTrees1,1);

j=1;
empty_val_index=1;
empty_val=[];
for i=1:num_within_seq1
    bol=isempty(all_qCompTrees1{i});
    if bol==1
        empty_val(empty_val_index)=i;
        empty_val_index=empty_val_index+1;
        continue;
    end
   aligned_q_tree_1{j}=all_qCompTrees1{i};
   j=j+1;  
end

num_within_seq1=numel(aligned_q_tree_1);

data_path2 = 'utils_data/NeuroData/maize02_seq/';
[all_qCompTrees2, all_compTrees2] = maize_load_from_mat_file(data_path2);
num_within_seq2=numel(all_qCompTrees2);
maize_show_all_shapes_sp_reg(all_compTrees2,2);

data_path3 = 'utils_data/NeuroData/maize03_seq_alt/';
[all_qCompTrees3, all_compTrees3] = maize_load_from_mat_file(data_path3);
num_within_seq3=numel(all_qCompTrees3);
maize_show_all_shapes_sp_reg(all_compTrees3,3);

%remove the null nodes
j=1;
empty_val_index=1;
empty_val=[];
for i=1:num_within_seq3
    bol=isempty(all_qCompTrees3{i});
    if bol==1
        empty_val(empty_val_index)=i;
        empty_val_index=empty_val_index+1;
        continue;
    end
   aligned_q_tree_3{j}=all_qCompTrees3{i};
   j=j+1;  
end

num_within_seq3=numel(aligned_q_tree_3);

data_path4 = 'utils_data/NeuroData/maize04_seq/';
[all_qCompTrees4, all_compTrees4] = maize_load_from_mat_file(data_path4);
num_within_seq4=numel(all_qCompTrees4);
maize_show_all_shapes_sp_reg(all_compTrees4,4);

data_path5 = 'utils_data/NeuroData/maize05_seq_alt/';
[all_qCompTrees5, all_compTrees5] = maize_load_from_mat_file(data_path5);
num_within_seq5=numel(all_qCompTrees5);
maize_show_all_shapes_sp_reg(all_compTrees5,5);

%remove the null nodes
j=1;
empty_val_index=1;
empty_val=[];
for i=1:num_within_seq5
    bol=isempty(all_qCompTrees5{i});
    if bol==1
        empty_val(empty_val_index)=i;
        empty_val_index=empty_val_index+1;
        continue;
    end
   aligned_q_tree_5{j}=all_qCompTrees5{i};
   j=j+1;  
end

num_within_seq5=numel(aligned_q_tree_5);

data_path6 = 'utils_data/NeuroData/maize06_seq/';
[all_qCompTrees6, all_compTrees6] = maize_load_from_mat_file(data_path6);
num_within_seq6=numel(all_qCompTrees6);
maize_show_all_shapes_sp_reg(all_compTrees6,6);

data_path7 = 'utils_data/NeuroData/maize07_seq_alt/';
[all_qCompTrees7, all_compTrees7] = maize_load_from_mat_file(data_path7);
num_within_seq7=numel(all_qCompTrees7);
maize_show_all_shapes_sp_reg(all_compTrees7,7);

%remove the null nodes
j=1;
empty_val_index=1;
empty_val=[];
for i=1:num_within_seq7
    bol=isempty(all_qCompTrees7{i});
    if bol==1
        empty_val(empty_val_index)=i;
        empty_val_index=empty_val_index+1;
        continue;
    end
   aligned_q_tree_7{j}=all_qCompTrees7{i};
   j=j+1;  
end

num_within_seq7=numel(aligned_q_tree_7);


aligned_q_tree_1{1}=aligned_q_tree_1{1};
%registration within_sequence
for i=2:num_within_seq1
    Q1=aligned_q_tree_1{i};
    Q2 = aligned_q_tree_1{i-1};
    [G,Q1p, Q2p] = ReparamPerm_qCompTrees_rad_4layers_v2(Q1, Q2, lam_m, lam_s, lam_p);
    aligned_q_tree_1{i}=Q1p;
    %aligned_q_tree_1{i-1}=Q2p;
end


aligned_q_tree_2{1}=all_qCompTrees2{1};
%registration within_sequence
for i=2:num_within_seq2
    Q1=all_qCompTrees2{i};
    Q2 = aligned_q_tree_2{i-1};
    [G,Q1p, Q2p] = ReparamPerm_qCompTrees_rad_4layers_v2(Q1, Q2, lam_m, lam_s, lam_p);
    aligned_q_tree_2{i}=Q1p;
    %aligned_q_tree_1{i-1}=Q2p;
end

aligned_q_tree_3{1}=aligned_q_tree_3{1};
%registration within_sequence
for i=2:num_within_seq3
    Q1=aligned_q_tree_3{i};
    Q2 = aligned_q_tree_3{i-1};
    [G,Q1p, Q2p] = ReparamPerm_qCompTrees_rad_4layers_v2(Q1, Q2, lam_m, lam_s, lam_p);
    aligned_q_tree_3{i}=Q1p;
    %aligned_q_tree_1{i-1}=Q2p;
end

aligned_q_tree_4{1}=all_qCompTrees4{1};
%registration within_sequence
for i=2:num_within_seq4
    Q1=all_qCompTrees4{i};
    Q2 = aligned_q_tree_4{i-1};
    [G,Q1p, Q2p] = ReparamPerm_qCompTrees_rad_4layers_v2(Q1, Q2, lam_m, lam_s, lam_p);
    aligned_q_tree_4{i}=Q1p;
    %aligned_q_tree_1{i-1}=Q2p;
end
aligned_q_tree_5{1}=aligned_q_tree_5{1};
%registration within_sequence
for i=2:num_within_seq5
    Q1=aligned_q_tree_5{i};
    Q2 = aligned_q_tree_5{i-1};
    [G,Q1p, Q2p] = ReparamPerm_qCompTrees_rad_4layers_v2(Q1, Q2, lam_m, lam_s, lam_p);
    aligned_q_tree_5{i}=Q1p;
    %aligned_q_tree_1{i-1}=Q2p;
end

aligned_q_tree_6{1}=all_qCompTrees6{1};
%registration within_sequence
for i=2:num_within_seq6
    Q1=all_qCompTrees6{i};
    Q2 = aligned_q_tree_6{i-1};
    [G,Q1p, Q2p] = ReparamPerm_qCompTrees_rad_4layers_v2(Q1, Q2, lam_m, lam_s, lam_p);
    aligned_q_tree_6{i}=Q1p;
    %aligned_q_tree_1{i-1}=Q2p;
end

aligned_q_tree_7{1}=aligned_q_tree_7{1};
%registration within_sequence
for i=2:num_within_seq7
    Q1=aligned_q_tree_7{i};
    Q2 = aligned_q_tree_7{i-1};
    [G,Q1p, Q2p] = ReparamPerm_qCompTrees_rad_4layers_v2(Q1, Q2, lam_m, lam_s, lam_p);
    aligned_q_tree_7{i}=Q1p;
    %aligned_q_tree_1{i-1}=Q2p;
end

[Mu,eigenVector,EVals,tree_pca,used_qCompTrees_Ready]=transfer_to_PCA_space_raw(aligned_q_tree_1,aligned_q_tree_2,aligned_q_tree_3,aligned_q_tree_4,aligned_q_tree_5,aligned_q_tree_6,aligned_q_tree_7,num_within_seq1,num_within_seq2,num_within_seq3,num_within_seq4,num_within_seq5,num_within_seq6,num_within_seq7);

%index of all raw trees
index=[1,4,5,10,11,14,15,19,20,23,24,30,31,35];

for i=1:2:14
    X=tree_pca(index(i):index(i+1),:)';
    X=ReSampleCurve(X,7);
      
    for k=1:7
       tree1=X(:,k)'*eigenVector';

       tree=tree1+Mu';

       recon{k}=tree;
    end

   for j=1:size(recon,2)
   
   if i==1
       aligned_q_tree{j}=unflattenCompTree_4layers_rad(recon{j}, lam_m, lam_s, lam_p, used_qCompTrees_Ready{1});
       aligned_q_tree_11{j}=aligned_q_tree{j};
       
   elseif i==3
       aligned_q_tree{j}=unflattenCompTree_4layers_rad(recon{j}, lam_m, lam_s, lam_p, used_qCompTrees_Ready{1});
       aligned_q_tree_22{j}=aligned_q_tree{j};
   elseif i==5
       aligned_q_tree{j}=unflattenCompTree_4layers_rad(recon{j}, lam_m, lam_s, lam_p, used_qCompTrees_Ready{1});
       aligned_q_tree_33{j}=aligned_q_tree{j};
   elseif i==7
       aligned_q_tree{j}=unflattenCompTree_4layers_rad(recon{j}, lam_m, lam_s, lam_p, used_qCompTrees_Ready{1});
       aligned_q_tree_44{j}=aligned_q_tree{j};
   elseif i==9
       aligned_q_tree{j}=unflattenCompTree_4layers_rad(recon{j}, lam_m, lam_s, lam_p, used_qCompTrees_Ready{1});
       aligned_q_tree_55{j}=aligned_q_tree{j};
   elseif i==11
       aligned_q_tree{j}=unflattenCompTree_4layers_rad(recon{j}, lam_m, lam_s, lam_p, used_qCompTrees_Ready{1});
       aligned_q_tree_66{j}=aligned_q_tree{j};
   elseif i==13
       aligned_q_tree{j}=unflattenCompTree_4layers_rad(recon{j}, lam_m, lam_s, lam_p, used_qCompTrees_Ready{1});
       aligned_q_tree_77{j}=aligned_q_tree{j};
   end
 end
end

%Across registration
start=1;
trees=[aligned_q_tree_11,aligned_q_tree_22,aligned_q_tree_33,aligned_q_tree_44,aligned_q_tree_55,aligned_q_tree_66,aligned_q_tree_77];

for i=1:7
    j=i+7;
    for k=1:6
    Q1=trees{j};
    Q2 = trees{i};
    [G,Q1p, Q2p] = ReparamPerm_qCompTrees_rad_4layers_v2(Q1, Q2, lam_m, lam_s, lam_p);
    trees{j}=Q1p;
    j=j+7;
    end
end

for i=1:7
    PC1{i}=qCompTree_to_CompTree_rad_4layers(trees{i});
    PC2{i}=qCompTree_to_CompTree_rad_4layers(trees{i+7});
    PC3{i}=qCompTree_to_CompTree_rad_4layers(trees{i+14});
    PC4{i}=qCompTree_to_CompTree_rad_4layers(trees{i+21});
    PC5{i}=qCompTree_to_CompTree_rad_4layers(trees{i+28});
    PC6{i}=qCompTree_to_CompTree_rad_4layers(trees{i+35});
    PC7{i}=qCompTree_to_CompTree_rad_4layers(trees{i+42});
end

gcf = figure;%('position', [0, 0, 1500, 800]);

set(gcf,'visible', 'on');
set(gcf, 'color', 'w');
view(0, 0);
axis equal; hold on;
box on; hold on;

set(gca,'XColor', 'none','YColor','none','ZColor','none');

maize_show_all_shapes_sp_reg(PC1,1);
maize_show_all_shapes_sp_reg(PC2,2);
maize_show_all_shapes_sp_reg(PC3,3);
maize_show_all_shapes_sp_reg(PC4,4);
maize_show_all_shapes_sp_reg(PC5,5);
maize_show_all_shapes_sp_reg(PC6,6);
maize_show_all_shapes_sp_reg(PC7,7);

%transfer to pca space
tNum=numel(trees);
[used_qCompTrees_Ready] = CompatMultiMax_rad_4layers(trees);

qX = [];

for i = 1:tNum
   qX(i, :) = flattenQCompTree_4layers_rad(used_qCompTrees_Ready{i}, lam_m, lam_s, lam_p);
end

qXT= qX';

[Mu,eigenVector, EVals] = performEigenAnalysis(qXT);


for i=1:tNum

sample{i}=qXT(:,i)-Mu;

tree_pca_curve(i,:)=sample{i}'* eigenVector;

end

%% compute modes of variation and random samples
%reshaping curves
start=1;
for i=1:7
    ending=start+6;
    X_curves(:,i)=reshape(tree_pca_curve(start:ending,:)',[],1);
    start=ending+1;
end


%PCA on curves
tic
[Mu_s, eigenVectors_s, EVals_s] = performEigenAnalysis(X_curves);
toc

alpha_1 = sqrt(EVals_s(1));
alpha_2 = sqrt(EVals_s(2));
alpha_3 = sqrt(EVals_s(3));
alpha_4 = sqrt(EVals_s(4));

Range=1.5; Step=0.5;

alpha_1_vec = -Range*alpha_2: (Step*alpha_2): Range*alpha_2;

mode_length= numel(alpha_1_vec);

tic
for i=1:mode_length
    qX_PC1{i} = alpha_1_vec(i)*eigenVectors_s(:,2) + Mu_s;
end
toc

X=zeros(7,48,7);

for i=1:mode_length
X_modes(i,:,:)=reshape(qX_PC1{i},[48 7]);
end


digit=(rand(1, 40))*1;

%for md=1:7
for md=1:5 %for rand sample
    %tic
    rand_sampleX =digit(md)*alpha_1 * eigenVectors_s(:,1) + digit(md)*alpha_2 *eigenVectors_s(:,2) + digit(md)*alpha_3 *eigenVectors_s(:,3) + Mu_s;
    %toc
    X1=reshape(rand_sampleX,[48 7]);%mean or random sampl
for i=1:size(X1,2) % for modes (X,3)
       %tree1=X_modes(md,:,i)*eigenVector'; % for modes
       tree1=X1(:,i)'*eigenVector'; 
       tree=tree1+Mu';

       recon{i}=tree;


end

for j=1:size(recon,2)

   PC1{j}=unflattenCompTree_4layers_rad(recon{j}, lam_m, lam_s, lam_p, used_qCompTrees_Ready{1});
   PC{j}=qCompTree_to_CompTree_rad_4layers(PC1{j});

end
%show_all_shapes_modes(PC,1);
maize_show_all_shapes_modes(PC,md);

end

return;
