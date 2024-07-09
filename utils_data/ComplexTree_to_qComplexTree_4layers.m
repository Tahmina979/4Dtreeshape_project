function qST = ComplexTree_to_qComplexTree_4layers(ST)

qST = SimpleTree_to_qSimpleTree(ST);

qST.q_children = cell(1, numel(ST.beta_children));
for i=1: numel(ST.beta_children)
    qST.q_children{i} = ComplexTree_to_qComplexTree(qST.beta_children{i});
%     qST.q_children{i}.q_children = qST.beta_children{i}.q_children;
end
    
end