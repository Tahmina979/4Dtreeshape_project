function [all_qCompTrees, all_compTrees] = load_neuroTrees_rad(data_path)

% --- get the swc files in the data_path ---

swc_files = get_filenames(data_path);
tree_num = numel(swc_files)

%% read data and extract apical dendrite
raw_trees = cell(1,tree_num);
ST = cell(1,tree_num);
qST = cell(1,tree_num);
% figure;
% axis equal; hold on;
all_compTrees = cell(1, tree_num);
for i=1:tree_num
    
    raw_trees{i} = read_swcdata( strcat( data_path, swc_files{i}) );
   
    all_compTrees{i} = compTree_from_swcdata_rad(raw_trees{i}, 4);  
    
end
%save("03_13.mat","all_compTrees");

all_qCompTrees = cell(1, tree_num);
for i= 1:tree_num
    all_qCompTrees{i} = CompTree_to_qCompTree_rad_4layers(all_compTrees{i});
end

end

