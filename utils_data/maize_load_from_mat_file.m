function [all_qCompTrees, all_compTrees] = maize_load_from_mat_file(data_path)


mat_files = dir(fullfile(data_path, '*.mat'));
tree_num = length(mat_files)

ST = cell(1,tree_num);
qST = cell(1,tree_num);
all_compTrees=cell(1,tree_num);
for i=1:tree_num
    tree=load(strcat(data_path,mat_files(i).name));
    all_compTrees{i} = tree.all_compTrees;          %for maize plants
end

all_qCompTrees = cell(1, tree_num);

for i= 1:tree_num
    
    bol=isempty(all_compTrees{i});
    if bol~=1
        all_qCompTrees{i} = CompTree_to_qCompTree_rad_4layers(all_compTrees{i});
    else
        all_qCompTrees{i}=[];

end

end

