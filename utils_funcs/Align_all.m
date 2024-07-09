function [aligned_trees]=Align_all (all_qCompTrees1,all_qCompTrees2,all_qCompTrees3,all_qCompTrees4,all_qCompTrees5,all_qCompTrees6,all_qCompTrees7)

num_within_seq1=numel(all_qCompTrees1);
num_within_seq2=numel(all_qCompTrees2);
num_within_seq3=numel(all_qCompTrees3);
num_within_seq4=numel(all_qCompTrees4);
num_within_seq5=numel(all_qCompTrees5);
num_within_seq6=numel(all_qCompTrees6);
num_within_seq7=numel(all_qCompTrees7);


aligned_q_tree_1=cell(1,num_within_seq1);
aligned_q_tree_2=cell(1,num_within_seq2);
aligned_q_tree_3=cell(1,num_within_seq3);
aligned_q_tree_4=cell(1,num_within_seq4);
aligned_q_tree_5=cell(1,num_within_seq5);
aligned_q_tree_6=cell(1,num_within_seq6);
aligned_q_tree_7=cell(1,num_within_seq7);

aligned_right=cell(1,7);


% parameters
lam_m = 1; 
lam_s = 1;
lam_p = 1;

aligned_q_tree_1{1}=all_qCompTrees1{1};
%registration within_sequence
for i=2:num_within_seq1
    Q1=all_qCompTrees1{i};
    Q2 = aligned_q_tree_1{i-1};
    [G,Q1p, Q2p] = ReparamPerm_qCompTrees_rad_4layers_v2(Q1, Q2, lam_m, lam_s, lam_p);
    aligned_q_tree_1{i}=Q1p;
    aligned_q_tree_1{i-1}=Q2p;
end

aligned_q_tree_2{1}=all_qCompTrees2{1};
for i=2:num_within_seq2
    Q1=all_qCompTrees2{i};
    Q2 = aligned_q_tree_2{i-1};
    [G,Q1p, Q2p] = ReparamPerm_qCompTrees_rad_4layers_v2(Q1, Q2, lam_m, lam_s, lam_p);
    aligned_q_tree_2{i}=Q1p;
    aligned_q_tree_2{i-1}=Q2p;
end

aligned_q_tree_3{1}=all_qCompTrees3{1};
for i=2:num_within_seq3
    Q1=all_qCompTrees3{i};
    Q2 = aligned_q_tree_3{i-1};
    [G,Q1p, Q2p] = ReparamPerm_qCompTrees_rad_4layers_v2(Q1, Q2, lam_m, lam_s, lam_p);
    aligned_q_tree_3{i}=Q1p;
    aligned_q_tree_3{i-1}=Q2p;
end

aligned_q_tree_4{1}=all_qCompTrees4{1};
for i=2:num_within_seq4
    Q1=all_qCompTrees4{i};
    Q2 = aligned_q_tree_4{i-1};
    [G,Q1p, Q2p] = ReparamPerm_qCompTrees_rad_4layers_v2(Q1, Q2, lam_m, lam_s, lam_p);
    aligned_q_tree_4{i}=Q1p;
    aligned_q_tree_4{i-1}=Q2p;
end

aligned_q_tree_5{1}=all_qCompTrees5{1};
for i=2:num_within_seq5
    Q1=all_qCompTrees5{i};
    Q2 = aligned_q_tree_5{i-1};
    [G,Q1p, Q2p] = ReparamPerm_qCompTrees_rad_4layers_v2(Q1, Q2, lam_m, lam_s, lam_p);
    aligned_q_tree_5{i}=Q1p;
    aligned_q_tree_5{i-1}=Q2p;
end

aligned_q_tree_6{1}=all_qCompTrees6{1};
for i=2:num_within_seq6
    Q1=all_qCompTrees6{i};
    Q2 = aligned_q_tree_6{i-1};
    [G,Q1p, Q2p] = ReparamPerm_qCompTrees_rad_4layers_v2(Q1, Q2, lam_m, lam_s, lam_p);
    aligned_q_tree_6{i}=Q1p;
    aligned_q_tree_6{i-1}=Q2p;
end

aligned_q_tree_7{1}=all_qCompTrees7{1};
for i=2:num_within_seq7
    Q1=all_qCompTrees7{i};
    Q2 = aligned_q_tree_7{i-1};
    [G,Q1p, Q2p] = ReparamPerm_qCompTrees_rad_4layers_v2(Q1, Q2, lam_m, lam_s, lam_p);
    aligned_q_tree_7{i}=Q1p;
    aligned_q_tree_7{i-1}=Q2p;
end


aligned_right{1}=aligned_q_tree_1{num_within_seq1};
aligned_right{2}=aligned_q_tree_2{num_within_seq2};
aligned_right{3}=aligned_q_tree_3{num_within_seq3};
aligned_right{4}=aligned_q_tree_4{num_within_seq4};
aligned_right{5}=aligned_q_tree_5{num_within_seq5};
aligned_right{6}=aligned_q_tree_6{num_within_seq6};
aligned_right{7}=aligned_q_tree_7{num_within_seq7};


%registration between sequence

for i=2:7

Q1=aligned_right{i};  
Q2=aligned_right{i-1};

[G,Q1p, Q2p] = ReparamPerm_qCompTrees_rad_4layers_v2(Q1, Q2, lam_m, lam_s, lam_p);
aligned_right{i}=Q1p;
aligned_right{i-1}=Q2p;

end

for i=7:2
Q2=aligned_right{i}; 
Q1=aligned_right{i-1};
[G,Q1p, Q2] = ReparamPerm_qCompTrees_rad_4layers_v2_within(Q1, Q2, lam_m, lam_s, lam_p);
aligned_right{i-1}=Q1p;
end

aligned_q_tree_1{num_within_seq1}=aligned_right{1}
aligned_q_tree_2{num_within_seq2}=aligned_right{2}
aligned_q_tree_3{num_within_seq3}=aligned_right{3}
aligned_q_tree_4{num_within_seq4}=aligned_right{4}
aligned_q_tree_5{num_within_seq5}=aligned_right{5}
aligned_q_tree_6{num_within_seq6}=aligned_right{6}
aligned_q_tree_7{num_within_seq7}=aligned_right{7}

%regiistration withiin sequences_1

for i=num_within_seq1:-1:2
    Q1=aligned_q_tree_1{i-1};
    Q2 = aligned_q_tree_1{i};
    [G,Q1p, Q2p] = ReparamPerm_qCompTrees_rad_4layers_v2_within(Q1, Q2, lam_m, lam_s, lam_p);
    aligned_q_tree_1{i-1}=Q1p;
    
    %aligned_q_tree_1{i}=Q2p;
end


%regiistration withiin sequences_2


for i=num_within_seq2:-1:2
    Q1=aligned_q_tree_2{i-1};
    Q2 = aligned_q_tree_2{i};
    [G,Q1p, Q2p] = ReparamPerm_qCompTrees_rad_4layers_v2_within(Q1, Q2, lam_m, lam_s, lam_p);
    aligned_q_tree_2{i-1}=Q1p;
   
    %aligned_q_tree_2{i}=Q2p;
end

%regiistration withiin sequences_3

for i=num_within_seq3:-1:2
    Q1=aligned_q_tree_3{i-1};
    Q2 = aligned_q_tree_3{i};
    [G,Q1p, Q2p] = ReparamPerm_qCompTrees_rad_4layers_v2_within(Q1, Q2, lam_m, lam_s, lam_p);
    aligned_q_tree_3{i-1}=Q1p;
    
    %aligned_q_tree_1{i}=Q2p;
end

%regiistration withiin sequences_4


for i=num_within_seq4:-1:2
    Q1=aligned_q_tree_4{i-1};
    Q2 = aligned_q_tree_4{i};
    [G,Q1p, Q2p] = ReparamPerm_qCompTrees_rad_4layers_v2_within(Q1, Q2, lam_m, lam_s, lam_p);
    aligned_q_tree_4{i-1}=Q1p;
   
    %aligned_q_tree_2{i}=Q2p;
end

%regiistration withiin sequences_5

for i=num_within_seq5:-1:2
    Q1=aligned_q_tree_5{i-1};
    Q2 = aligned_q_tree_5{i};
    [G,Q1p, Q2p] = ReparamPerm_qCompTrees_rad_4layers_v2_within(Q1, Q2, lam_m, lam_s, lam_p);
    aligned_q_tree_5{i-1}=Q1p;
    
    %aligned_q_tree_1{i}=Q2p;
end



%regiistration withiin sequences_6


for i=num_within_seq6:-1:2
    Q1=aligned_q_tree_6{i-1};
    Q2 = aligned_q_tree_6{i};
    [G,Q1p, Q2p] = ReparamPerm_qCompTrees_rad_4layers_v2_within(Q1, Q2, lam_m, lam_s, lam_p);
    aligned_q_tree_6{i-1}=Q1p;
   
    %aligned_q_tree_2{i}=Q2p;
end

%regiistration withiin sequences_7


for i=num_within_seq7:-1:2
    Q1=aligned_q_tree_7{i-1};
    Q2 = aligned_q_tree_7{i};
    [G,Q1p, Q2p] = ReparamPerm_qCompTrees_rad_4layers_v2_within(Q1, Q2, lam_m, lam_s, lam_p);
    aligned_q_tree_7{i-1}=Q1p;
   
    %aligned_q_tree_2{i}=Q2p;
end


aligned_trees=[aligned_q_tree_1,aligned_q_tree_2,aligned_q_tree_3,aligned_q_tree_4,aligned_q_tree_5,aligned_q_tree_6,aligned_q_tree_7];
save('aligned_trees.mat',"aligned_trees");
end