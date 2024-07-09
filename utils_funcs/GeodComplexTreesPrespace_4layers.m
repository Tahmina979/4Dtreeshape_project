function [A, qA] = GeodComplexTreesPrespace_4layers(Q1, Q2, geo_param, start_origin)


[~,qA] = GeodSimpleTreesPrespace(Q1, Q2, geo_param);
for i=1: numel(qA)
    qA{i}.q_children = cell(1, numel(Q1.q_children));
end


for i=1: numel(Q1.q_children)
    [~,qA_children] = GeodComplexTreesPrespace(Q1.q_children{i}, Q2.q_children{i}, geo_param);
    for j = 1: numel(qA_children)
        qA{j}.q_children{i} = qA_children{j};
    end
end



% convert to AC-space (original trees)
A = cell(1, numel(qA));
stp = geo_param;
for i=1:numel(qA)
    A{i} = qComplexTree_to_ComplexTree_4layers(qA{i});
end
