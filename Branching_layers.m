clc;
clear;
close all;

addpath('utils_data','utils_draw');

gcf = figure;
set(gcf,'visible', 'on');
set(gcf, 'color', 'w');
view(0, 0);
axis equal; hold on;
box on; hold on;

set(gca,'XColor', 'none','YColor','none','ZColor','none');

data_path1= 'utils_data/NeuroData/maize01_seq/';
[all_qCompTrees1, all_compTrees1] = load_from_mat_file(data_path1);
show_all_shapes_reg(all_compTrees1,1);
