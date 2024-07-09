function [tree1, tree2]=Align (all_qCompTrees1,data_path1,all_qCompTrees2,data_path2)

num_within_seq1=numel(all_qCompTrees1);
num_within_seq2=numel(all_qCompTrees2);


aligned_q_tree_1=cell(1,num_within_seq1);
tree1=cell(1,num_within_seq1);

aligned_q_tree_2=cell(1,num_within_seq2);
tree2=cell(1,num_within_seq2);

% parameters
lam_m = 1; 
lam_s = 1;
lam_p = 1;

%registration between sequence
Q1=all_qCompTrees2{1};
Q2=all_qCompTrees1{1};

[G,Q1p, Q2p] = ReparamPerm_qCompTrees_rad_4layers_v2(Q1, Q2, lam_m, lam_s, lam_p);
aligned_q_tree_2{1}=Q1p;
aligned_q_tree_1{1}=Q2p;


for i=2:num_within_seq1
    Q1=all_qCompTrees1{i};
    Q2 = aligned_q_tree_1{i-1};
    [G,Q1p, Q2p] = ReparamPerm_qCompTrees_rad_4layers_v2(Q1, Q2, lam_m, lam_s, lam_p);
    aligned_q_tree_1{i}=Q1p;
    aligned_q_tree_1{i-1}=Q2p;
end


for i=2:num_within_seq2
    Q1=all_qCompTrees2{i};
    Q2 = aligned_q_tree_2{i-1};
    [G,Q1p, Q2p] = ReparamPerm_qCompTrees_rad_4layers_v2(Q1, Q2, lam_m, lam_s, lam_p);
    aligned_q_tree_2{i}=Q1p;
    aligned_q_tree_2{i-1}=Q2p;
end


%registration between sequence
Q1=aligned_q_tree_2{num_within_seq2};
Q2=aligned_q_tree_1{num_within_seq1};

[G,Q1p, Q2p] = ReparamPerm_qCompTrees_rad_4layers_v2(Q1, Q2, lam_m, lam_s, lam_p);
aligned_q_tree_2{num_within_seq2}=Q1p;
aligned_q_tree_1{num_within_seq1}=Q2p;

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
tree1=aligned_q_tree_1;
tree2=aligned_q_tree_2;

end