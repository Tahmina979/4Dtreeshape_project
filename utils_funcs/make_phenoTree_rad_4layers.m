function Q = make_phenoTree_rad_4layers(beta0,beta_children,tk_sideLocs,b00, beta0_rad, beta_rad)

if nargin<4
    Q.b00_startP = [0;0;0;0];
else
    Q.b00_startP = b00;
end

Q.beta0 = beta0;
Q.beta0_rad = beta0_rad;
Q.beta_children = beta_children;
Q.beta_rad = beta_rad;

Q.K_sideNum = numel(beta_children);
[Q.dimension,Q.T0_pointNum] = size(beta0);

Q.t_paras = linspace(0,1,Q.T0_pointNum);

Q.len0 = Q.t_paras(end);


Q.tk_sideLocs = tk_sideLocs;

Q.T_sidePointNums = zeros(1, Q.K_sideNum);
Q.len = zeros(1, Q.K_sideNum);
for k=1:Q.K_sideNum
    Q.beta{k} = Q.beta_children{k}.beta0;
    Q.T_sidePointNums(k) = size(Q.beta{k}, 2);
    t = linspace(0, 1, Q.T_sidePointNums(k));
    Q.len(k) = trapz( t, sum(Q.beta{k}.^2,1), 2);
end

Q = orderfields(Q,{'t_paras', 'beta0','T0_pointNum','len0', 'K_sideNum','beta','T_sidePointNums','len', ...
    'tk_sideLocs','dimension','b00_startP', 'beta0_rad', 'beta_rad','beta_children'});