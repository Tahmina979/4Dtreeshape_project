function ST = qSimpleTree_to_SimpleTree(qST)

ST = qST;
% ST = rmfield(qST, {'q','q0','b00','len','len0','s','sk'});

ST.beta0 = q_to_curve(qST.q0) + repmat(qST.b00_startP,1,qST.T0_pointNum);

if isempty(ST.tk_sideLocs) == 0
    betak0 = interp1(ST.t_paras, ST.beta0', ST.tk_sideLocs)';
else
    betak0 = [];
end

ST.beta = cell(1,ST.K_sideNum);
for k=1:ST.K_sideNum
    
    if isempty(qST.q{k})
        ST.beta{k} = [0, 0, 0]'+ + repmat(betak0(:,k), 1, size(qST.q{k}, 2));
    else
        ST.beta{k} = q_to_curve(qST.q{k});   % --- now the starting point of ST.beta{k} is at the orgin.
        ST.beta{k} = ST.beta{k} + repmat(betak0(:,k), 1, size(qST.q{k}, 2));
%     end
    end


% ST = orderfields(ST, {'t','beta0','T0','K','beta','T','tk','d'});

end