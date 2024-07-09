function [A, qA] = GeodComplexTreesPrespace(Q1, Q2, geo_param, start_origin, weights)

if nargin < 5,
    weights = [1 1 1];
end

if nargin < 4,
    start_origin = true;
end

[~,qA] = GeodSimpleTreesPrespace(Q1, Q2, geo_param, start_origin, weights); 
for i=1: numel(qA)
    qA{i}.q_children = cell(1, numel(Q1.q_children));
end

% if isfield(Q1, 'q_children_p') == 0
%     Q1.q_children_p = cell(1, 0);
% end
% 
% if isfield(Q2, 'q_children_p') == 0
%     Q2.q_children_p = cell(1, 0);
% end


for i=1: numel(Q1.q_children)
    [~, qA_children] = GeodSimpleTreesPrespace(Q1.q_children{i}, Q2.q_children{i}, geo_param, start_origin, weights);
    for j = 1: numel(qA_children)
        qA{j}.q_children{i} = qA_children{j};
    end
end

% convert to AC-space (original trees)
A = cell(1, numel(qA));
stp = geo_param;
for i=1:numel(qA)
    A{i} = qComplexTree_to_ComplexTree(qA{i});
end

% --- rewritten by me ---
% A{stp} = qSimpleTree_to_SismpleTree(Q2);

