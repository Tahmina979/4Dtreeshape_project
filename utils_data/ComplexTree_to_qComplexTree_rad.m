function qST = ComplexTree_to_qComplexTree_rad(ST)

qST = SimpleTree_to_qSimpleTree_rad(ST);

qST.q_children = cell(1, numel(ST.beta_children));
for i=1: numel(ST.beta_children)
    qST.q_children{i} = SimpleTree_to_qSimpleTree_rad(qST.beta_children{i});
    
    qST.q_children{i}.q_children = cell(1, numel(qST.q_children{i}.q));
    
    for j=1: numel(qST.q_children{i}.q)
        qST.q_children{i}.q_children{j}.q0 = qST.q_children{i}.q{j};
        
%         qST.q_children{i}.q_children{j}.t_paras = Compute_t_paras(qST.q_children{i}.q_children{j}.q0);
        qST.q_children{i}.q_children{j}.func_t_parasAndRad ......
                            = qST.beta_children{i}.beta_children{j}.func_t_parasAndRad;
                        
        qST.q_children{i}.q_children{j}.data_t_parasAndRad ......
            = qST.beta_children{i}.beta_children{j}.data_t_parasAndRad;
        
        qST.q_children{i}.q_children{j}.t_paras ......
            = qST.beta_children{i}.beta_children{j}.t_paras;
    end
end
    
end