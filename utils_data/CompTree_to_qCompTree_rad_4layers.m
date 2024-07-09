function qCompTree_4Ls = CompTree_to_qCompTree_rad_4layers(compTree_4Ls)

qCompTree_4Ls = ST_to_qST_rad(compTree_4Ls);

qCompTree_4Ls.q_children = cell(1, numel(compTree_4Ls.beta_children));
for i=1: numel(compTree_4Ls.beta_children)
    qCompTree_4Ls.q_children{i} = CompTree_to_qCompTree_rad_3layers(compTree_4Ls.beta_children{i});
%     qST.q_children{i}.q_children = qST.beta_children{i}.q_children;
end
    
end