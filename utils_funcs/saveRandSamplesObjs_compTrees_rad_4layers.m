% clc
% close all;
% gcf = figure('position', [100, 100, 1000, 600]);
% set(gcf,'visible', 'on');
% set(gcf, 'color', 'w');
% view(0, 90);
% axis equal; hold on;

% A10 is the geodesic
function [folder_name] = saveRandSamplesObjs_compTrees_rad_4layers(A10, data_path, idx1, idx2)

CT = A10;

addpath Yanirk_code

% c = distinguishable_colors(50,[.9,.9,.9;0,0,0],[20,20]);

c = distinguishable_colors(70);

used_idxes = -1;
folder_name = ['randSamplesOBJs_', data_path];
mkdir(folder_name);

for i =1: numel(CT)
    
%     y_added = ceil(i/19) * -0;
%     x_added = (mod(i, 19)+1) * 300;
%     
%     clear first_pt_x first_pt_y first_pt_z
%     
%     first_pt_x = (CT{i}.beta0(1, 1));
%     first_pt_y = (CT{i}.beta0(2, 1));
%     first_pt_z = (CT{i}.beta0(3, 1));
    
    compTree2obj_rad_4layer(CT{i}, [folder_name,'/',num2str(i)]);
end
    