function qST = ST_to_qST_rad(ST)

qST = rmfield(ST,{'beta0','beta'});

if isfield(qST, 'beta_children')
    qST = rmfield(qST,{'beta_children'});
end

qST.q0 = curve_to_q(ST.beta0);
qST.s = cumtrapz(qST.t_paras, sum(qST.q0.^2,1) );
qST.len0 = qST.s(end);
qST.s = qST.s/qST.len0;
qST.sk = interp1(qST.t_paras, qST.s, qST.tk_sideLocs);

qST.q = cell(1,ST.K_sideNum);
qST.len = zeros(1,ST.K_sideNum);
for k=1:qST.K_sideNum
    qST.q{k} = curve_to_q(ST.beta{k});
    qST.len(k) = trapz( sum(qST.q{k}.^2,1) / (qST.T_sidePointNums(k)-1));
end

qST.b00_startP = ST.beta0(:,1);

% if isfield(qST, 'beta_children')
%     qST = orderfields(qST, ...
%         {'t_paras','s','q0','T0_pointNum','len0','K_sideNum','q','T_sidePointNums','len','tk_sideLocs',...
%         'sk','dimension','b00_startP', 'beta_children', 'beta0_rad', 'beta_rad'});
% else
qST = orderfields(qST, ...
    {'t_paras','s','q0','T0_pointNum','len0','K_sideNum','q','T_sidePointNums','len','tk_sideLocs',...
    'sk','dimension','b00_startP', 'beta0_rad', 'beta_rad'});
% end

end