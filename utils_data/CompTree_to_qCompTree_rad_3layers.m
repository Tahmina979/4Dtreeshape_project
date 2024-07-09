function qCompTree_3Ls = CompTree_to_qCompTree_rad_3layers(compTree_3Ls)

qCompTree_3Ls = ST_to_qST_rad(compTree_3Ls);

qCompTree_3Ls.q_children = cell(1, numel(compTree_3Ls.beta_children));
for i=1: numel(compTree_3Ls.beta_children)
    qCompTree_3Ls.q_children{i} = ST_to_qST_rad(compTree_3Ls.beta_children{i});
end
    
end