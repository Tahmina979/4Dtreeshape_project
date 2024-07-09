function [ST_1] = trunkRep_to_compTree_4layers_rad(trunk1)

%   
ST_1 = trunkRep_to_ST_rad(trunk1);

ST_1.beta_children = cell(1, numel(trunk1.children));

for i = 1: numel(trunk1.children)
    ST_1.beta_children{i} = trunkRep_to_compTree_3layers_rad(trunk1.children{i});

end

end

