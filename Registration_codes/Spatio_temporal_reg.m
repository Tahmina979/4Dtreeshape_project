clc;
clear;
close all;

% This code will take time (arount 2-3 hour to register all the shapes),
% The time depends on the number of sequences we are trying to align and
% the totoal number of 3D trees in a sequence
addpath('utils_data','utils_draw',"OpenCurvesRn",'utils_funcs',"utils_statModels");

lam_m = 1; 
lam_s = 1;
lam_p = 1;

total_seq=7;
%for tomato plants, we have to use respective directory


ctree_1 = load('Registered_maize_plants/seq1_alligned.mat');

ctree_2 = load('Registered_maize_plants/seq2_alligned.mat');

ctree_3 = load('Registered_maize_plants/seq3_alligned.mat');

ctree_4 = load('Registered_maize_plants/seq4_alligned.mat');

ctree_5 = load('Registered_maize_plants/seq5_alligned.mat');

ctree_6 = load('Registered_maize_plants/seq6_alligned.mat');

ctree_7 = load('Registered_maize_plants/seq7_alligned.mat');

num_within_seq1=numel(ctree_1.aligned_q_tree_1);
num_within_seq2=numel(ctree_2.aligned_q_tree_2);
num_within_seq3=numel(ctree_3.aligned_q_tree_3);
num_within_seq4=numel(ctree_4.aligned_q_tree_4);
num_within_seq5=numel(ctree_5.aligned_q_tree_5);
num_within_seq6=numel(ctree_6.aligned_q_tree_6);
num_within_seq7=numel(ctree_7.aligned_q_tree_7);


[Mu,eigenVector,EVals,tree_pca,used_qCompTrees_Ready]=transfer_to_PCA_space(ctree_1,ctree_2,ctree_3,ctree_4,ctree_5,ctree_6,ctree_7,num_within_seq1,num_within_seq2,num_within_seq3,num_within_seq4,num_within_seq5,num_within_seq6,num_within_seq7);
tic
%X1=tree_pca(1:10,:)';
X1=tree_pca(1:num_within_seq1,:)'; %first sequence
X1=ReSampleCurve(X1,7); %seven is the highest length of all curves, so we interpolate all into ten, For tomato it is ten

X(:,1)=reshape(X1,[],1);

index=[num_within_seq1+1,num_within_seq1+num_within_seq2,num_within_seq1+num_within_seq2+1,num_within_seq1+num_within_seq2+num_within_seq3,num_within_seq1+num_within_seq2+num_within_seq3+1,num_within_seq1+num_within_seq2+num_within_seq3+num_within_seq4,num_within_seq1+num_within_seq2+num_within_seq3+num_within_seq4+1,num_within_seq1+num_within_seq2+num_within_seq3+num_within_seq4+num_within_seq5,num_within_seq1+num_within_seq2+num_within_seq3+num_within_seq4+num_within_seq5+1,num_within_seq1+num_within_seq2+num_within_seq3+num_within_seq4+num_within_seq5+num_within_seq6,num_within_seq1+num_within_seq2+num_within_seq3+num_within_seq4+num_within_seq5+num_within_seq6+1,num_within_seq1+num_within_seq2+num_within_seq3+num_within_seq4+num_within_seq5+num_within_seq6+num_within_seq7];
j=2;
for i=1:2:12
    X2=tree_pca(index(i):index(i+1),:)';
    X2=ReSampleCurve(X2,7); %seven is the highest length of all curves, so we interpolate all into ten, For tomato it is ten
    % temporal registration, it can be done later
    [dist,X2n,q2n,X1,X2,q1, gamI, distbefore]=mygeod(X1, X2);
    X(:,j)=reshape(X2n,[],1);
    j=j+1;
  
end

start=1;
ending=7; % for tomato 10
for i =1:total_seq 
    seq(:,start:ending)=reshape(X(:,i),[size(X1,1),size(X1,2)]);
    %start=start+10;
    %ending=ending+10;
    start=start+7;
    ending=ending+7;
end

start=1;
ending=size(X1,2);
for j=1:total_seq
    for i=start:ending
      
       tree1=seq(:,i)'*eigenVector';

       tree=tree1+Mu';

       recon{i}=tree;

       sequence{i}=unflattenCompTree_4layers_rad(recon{i}, lam_m, lam_s, lam_p, used_qCompTrees_Ready{1});
      
    end 
    start=start+7;
    ending=ending+7;
    %start=start+10;
    %ending=ending+10;
end
toc

tic

% Across spatial registration
count=0;
for i=1:7 %for tomato 10
    control=i;
    Q2=sequence{i};
    for j=1:total_seq-1
        %control=control+10; % for tomato
        control=control+7; % for maize
        Q1=sequence{control};
        [G,Q1p, Q2p] = ReparamPerm_qCompTrees_rad_4layers_v2(Q1, Q2, lam_m, lam_s, lam_p);
        sequence{control}=Q1p;
        count=count+1
    end
end

% for tomato it will be 10
num_within_seq1=7
num_within_seq2=7;
num_within_seq3=7;
num_within_seq4=7;
num_within_seq5=7;
num_within_seq6=7;
num_within_seq7=7;

toc
for i=1:num_within_seq1
  aligned_q_tree_1{i}=sequence{i};
end

p=num_within_seq1+1; p_end=num_within_seq1+num_within_seq2;

for i=p:p_end
  aligned_q_tree_2{i-p+1}=sequence{i};
end

p=p_end+1; p_end=p_end+num_within_seq3;
for i=p:p_end
   aligned_q_tree_3{i-p+1}=sequence{i};
end
p=p_end+1; p_end=p_end+num_within_seq4;
for i=p:p_end

    aligned_q_tree_4{i-p+1}=sequence{i};
end
p=p_end+1; p_end=p_end+num_within_seq5;
for i=p:p_end

   aligned_q_tree_5{i-p+1}=sequence{i};
end

p=p_end+1; p_end=p_end+num_within_seq6;
for i=p:p_end
 aligned_q_tree_6{i-p+1}=sequence{i};
end
p=p_end+1; p_end=p_end+num_within_seq7;
for i=p:p_end
   aligned_q_tree_7{i-p+1}=sequence{i};
end

save('seq1s_alligned','aligned_q_tree_1');
save('seq2s_alligned','aligned_q_tree_2');
save('seq3s_alligned','aligned_q_tree_3');
save('seq4s_alligned','aligned_q_tree_4');
save('seq5s_alligned','aligned_q_tree_5');
save('seq6s_alligned','aligned_q_tree_6');
save('seq7s_alligned','aligned_q_tree_7');

return;