clc;
clear;
close all;

data_path = 'utils_data/tomato_1/';
[all_qCompTrees, compTrees] = load_neuroTrees_rad(data_path);

% all_qCompTrees are the structered trees with layer information and we
% save those as .mat files


