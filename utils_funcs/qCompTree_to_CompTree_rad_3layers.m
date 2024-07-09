function ST = qCompTree_to_CompTree_rad_3layers(qST)

ST = qST_to_ST_rad(qST);
ST.beta_children = cell(1, numel(qST.q_children));

for i=1:numel(qST.q_children)
    if numel(ST.beta{i}) ~= 0
        qST.q_children{i}.b00_startP = ST.beta{i}(:, 1);
    else
       qST.q_children{i}.b00_startP = [0;0;0]; 
    end
    
    ST.beta_children{i} = qST_to_ST_rad(qST.q_children{i});
end



end