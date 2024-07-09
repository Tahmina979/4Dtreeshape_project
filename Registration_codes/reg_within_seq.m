clc;
clear;
addpath("utils_data", "utils_funcs", "utils_statModels", "utils_draw")

data_path1 = 'utils_data/NeuroData/maize01_seq/';
[all_qCompTrees1, all_compTrees1] = maize_load_from_mat_file(data_path1);
show_shapes(all_compTrees1);
return;
data_path2 = 'utils_data/NeuroData/maize02_seq/';
[all_qCompTrees2, all_compTrees2] = maize_load_from_mat_file(data_path2);
show_shapes(all_compTrees2);

data_path3 = 'utils_data/NeuroData/maize03_seq/';
[all_qCompTrees3, all_compTrees3] = maize_load_from_mat_file(data_path3);
show_shapes(all_compTrees3);
data_path4 = 'utils_data/NeuroData/maize04_seq/';
[all_qCompTrees4, all_compTrees4] = maize_load_from_mat_file(data_path4);
show_shapes(all_compTrees4);
data_path5 = 'utils_data/NeuroData/maize05_seq/';
[all_qCompTrees5, all_compTrees5] = maize_load_from_mat_file(data_path5);
show_shapes(all_compTrees5);
data_path6 = 'utils_data/NeuroData/maize06_seq/';
[all_qCompTrees6, all_compTrees6] = maize_load_from_mat_file(data_path6);
show_shapes(all_compTrees6);
data_path7 = 'utils_data/NeuroData/maize07_seq/';
[all_qCompTrees7, all_compTrees7] = maize_load_from_mat_file(data_path7);
show_shapes(all_compTrees7);


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
    %aligned_q_tree_1{i-1}=Q2p;
end

aligned_q_tree_2{1}=all_qCompTrees2{1};
for i=2:num_within_seq2
    Q1=all_qCompTrees2{i};
    Q2 = aligned_q_tree_2{i-1};
    [G,Q1p, Q2p] = ReparamPerm_qCompTrees_rad_4layers_v2(Q1, Q2, lam_m, lam_s, lam_p);
    aligned_q_tree_2{i}=Q1p;
   % aligned_q_tree_2{i-1}=Q2p;
end

aligned_q_tree_3{1}=all_qCompTrees3{1};
for i=2:num_within_seq3
    Q1=all_qCompTrees3{i};
    Q2 = aligned_q_tree_3{i-1};
    [G,Q1p, Q2p] = ReparamPerm_qCompTrees_rad_4layers_v2(Q1, Q2, lam_m, lam_s, lam_p);
    aligned_q_tree_3{i}=Q1p;
    %aligned_q_tree_3{i-1}=Q2p;
end

aligned_q_tree_4{1}=all_qCompTrees4{1};
for i=2:num_within_seq4
    Q1=all_qCompTrees4{i};
    Q2 = aligned_q_tree_4{i-1};
    [G,Q1p, Q2p] = ReparamPerm_qCompTrees_rad_4layers_v2(Q1, Q2, lam_m, lam_s, lam_p);
    aligned_q_tree_4{i}=Q1p;
   % aligned_q_tree_4{i-1}=Q2p;
end

aligned_q_tree_5{1}=all_qCompTrees5{1};
for i=2:num_within_seq5
    Q1=all_qCompTrees5{i};
    Q2 = aligned_q_tree_5{i-1};
    [G,Q1p, Q2p] = ReparamPerm_qCompTrees_rad_4layers_v2(Q1, Q2, lam_m, lam_s, lam_p);
    aligned_q_tree_5{i}=Q1p;
    %aligned_q_tree_5{i-1}=Q2p;
end

aligned_q_tree_6{1}=all_qCompTrees6{1};
for i=2:num_within_seq6
    Q1=all_qCompTrees6{i};
    Q2 = aligned_q_tree_6{i-1};
    [G,Q1p, Q2p] = ReparamPerm_qCompTrees_rad_4layers_v2(Q1, Q2, lam_m, lam_s, lam_p);
    aligned_q_tree_6{i}=Q1p;
    %aligned_q_tree_6{i-1}=Q2p;
end

aligned_q_tree_7{1}=all_qCompTrees7{1};
for i=2:num_within_seq7
    Q1=all_qCompTrees7{i};
    Q2 = aligned_q_tree_7{i-1};
    [G,Q1p, Q2p] = ReparamPerm_qCompTrees_rad_4layers_v2(Q1, Q2, lam_m, lam_s, lam_p);
    aligned_q_tree_7{i}=Q1p;
   % aligned_q_tree_7{i-1}=Q2p;
end

save("seq1_alligned.mat","aligned_q_tree_1");
save("seq2_alligned.mat","aligned_q_tree_2");
save("seq3_alligned.mat","aligned_q_tree_3");
save("seq4_alligned.mat","aligned_q_tree_4");
save("seq5_alligned.mat","aligned_q_tree_5");
save("seq6_alligned.mat","aligned_q_tree_6");
save("seq7_alligned.mat","aligned_q_tree_7");
