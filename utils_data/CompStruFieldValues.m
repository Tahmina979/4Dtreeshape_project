function B = CompStruFieldValues(B, beta, children_inds, id)

[d,T0] = size(beta{id});
t = linspace(0,1,T0);

K = numel(beta(children_inds{id}));

T = zeros(1,K);
tk = zeros(1,K);

for m=1:K
    T(m) = size(beta{m},2);
    bdist = sum( (beta{id}-repmat(beta{children_inds{id}(m)}(:,1),1,T0)).^2, 1 );
    [~,ii] = min(bdist);
    tk(m) = t(ii);
end
B.t_paras = t;
B.beta0 = beta{id};
B.T0_pointNum = T0;
B.K_sideNum = K;
B.beta = beta(children_inds{id});
B.T_sidePointNums = T;
B.tk_sideLocs = tk;
B.dimension = d;
B.beta_children = cell(1, K);

% --- refine the tk ---



end
