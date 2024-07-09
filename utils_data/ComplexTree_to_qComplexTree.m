function qST = ComplexTree_to_qComplexTree(ST)

qST = SimpleTree_to_qSimpleTree(ST);

qST.q_children = cell(1, numel(ST.beta_children));
for i=1: numel(ST.beta_children)
    qST.q_children{i} = SimpleTree_to_qSimpleTree(qST.beta_children{i});
end
    
end