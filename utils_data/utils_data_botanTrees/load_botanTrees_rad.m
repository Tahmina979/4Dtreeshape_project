function [all_qCompTrees, all_compTrees] = load_botanTrees_rad(folder_name)

addpath('utils_data');
F = dir([folder_name, '/*.txtskl']);
for i=1: size(F, 1)
    s =  F(i).name;
%     s(find(isspace(s))) = [];
    Files{i} = s;
end

tNum = size(F, 1);

for i_tree = 1: tNum
    
    [branch{i_tree}, branch_num{i_tree}] = ReadSkelFile([folder_name, '/', Files{i_tree}]);
        
    if numel(branch_num{i_tree}) == 1
        branch_num{i_tree}(2) = 0;
        branch{i_tree}(2, 1) = branch{i_tree}(1, 1);
    end
    
    if numel(branch_num{i_tree}) == 2
        branch_num{i_tree}(3) = 0;
        branch{i_tree}(3, 1) = branch{i_tree}(1, 1);
        
        branch_num{i_tree}(4) = 0;
        branch{i_tree}(4, 1) = branch{i_tree}(1, 1);
    end
    
    if numel(branch_num{i_tree}) == 3
        branch_num{i_tree}(4) = 0;
        branch{i_tree}(4, 1) = branch{i_tree}(1, 1);
    end
    
    clear trunk1 branch1
    trunk1 = TrunkStruction(branch{i_tree}, branch_num{i_tree});
    branch1 = branch{i_tree};
    
    % --- add the 4th layer to the trunk ---
    for i=1: numel(trunk1.children)
        for j=1: numel(trunk1.children{i}.children)
            trunk1.children{i}.children{j}.children = [];
            for k =1: numel(trunk1.children{i}.children{j}.bifurcation)
                Layer4_c_id = trunk1.children{i}.children{j}.bifurcation{k}.child;% �������ֻȡ��һ���� just for safe.
                trunk1.children{i}.children{j}.children{k} = branch1(4, Layer4_c_id);
    %             trunk1.children{i}.children{j}.children{k}.bifurcation = [];
            end
        end
    end

    % --- adjust the data structure for the 4-th layer
    for i=1: numel(trunk1.children)
        for ii = 1: numel(trunk1.children{i}.children)   
        k = 1;
        new_children = { };
        new_bifurcation = { };
        for j = 1: numel(trunk1.children{i}.children{ii}.children)
            if numel(trunk1.children{i}.children{ii}.children{j}) == 1
                new_children{k} = trunk1.children{i}.children{ii}.children{j};
                new_bifurcation{k} = trunk1.children{i}.children{ii}.bifurcation{j};
                k = k+1;
            elseif numel(trunk1.children{i}.children{ii}.children{j}) ==2
                new_children{k} = trunk1.children{i}.children{ii}.children{j}(1);
                new_bifurcation{k} = trunk1.children{i}.children{ii}.bifurcation{j};
                k = k+1;
                new_children{k} = trunk1.children{i}.children{ii}.children{j}(2);
                new_bifurcation{k} = trunk1.children{i}.children{ii}.bifurcation{j};
                k = k+1;
            end
        end

        trunk1.children{i}.children{ii}.children = new_children;
        trunk1.children{i}.children{ii}.bifurcation = new_bifurcation;

        end
    end

    % --- add the 'bifurcation' value to the 4-th layer
    for i=1: numel(trunk1.children)
        for j=1: numel(trunk1.children{i}.children)
            for k =1: numel(trunk1.children{i}.children{j}.children)
                trunk1.children{i}.children{j}.children{k}.bifurcation = [];
            end
        end
    end
    % --- adding ends here 

    % --- calculate the t_value for the third layer ---
    for i=1: numel(trunk1.children)
        for j=1: numel(trunk1.children{i}.children)
            trunk1.children{i}.children{j} = calcu_t_value_1_branch(trunk1.children{i}.children{j});
        end
    end
    
    trunk{i_tree} = trunk1;
     
end

all_compTrees= cell(1, tNum);
all_qCompTrees= cell(1, tNum);
for i= 1:tNum
    all_compTrees{i} = trunkRep_to_compTree_4layers_rad(trunk{i});
    all_qCompTrees{i} = CompTree_to_qCompTree_rad_4layers(all_compTrees{i});
end

end

