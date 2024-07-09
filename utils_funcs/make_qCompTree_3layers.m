function Q = make_qCompTree_3layers(q0,q_children,sk,b00)

if nargin<4
    Q.b00_startP = [0;0;0:0];
else
    Q.b00_startP = b00;
end

Q.q0 = q0;
Q.q_children = q_children;
Q.sk = sk;

Q.K_sideNum = numel(q_children);
[Q.dimension,Q.T0_pointNum] = size(q0);

Q.t_paras = linspace(0,1,Q.T0_pointNum);

Q.s = cumtrapz( Q.t_paras, sum(Q.q0.^2,1), 2);
Q.s = Q.s + linspace(0, 1e-4, length(Q.s)); % prevent Q.s to be zero vector
Q.len0 = Q.s(end);
Q.s = Q.s/Q.len0;

Q.tk_sideLocs = interp1(Q.s,Q.t_paras, Q.sk);

Q.T_sidePointNums = zeros(1, Q.K_sideNum);
Q.len = zeros(1, Q.K_sideNum);
Q.q = cell(1, 0);
for k=1:Q.K_sideNum
    Q.q{k} = Q.q_children{k}.q0;
    Q.T_sidePointNums(k) = size(Q.q{k}, 2);
    t = linspace(0, 1, Q.T_sidePointNums(k));
    Q.len(k) = trapz( t, sum(Q.q{k}.^2,1), 2);
end

Q = orderfields(Q,{'t_paras','s','q0','T0_pointNum','len0',...
    'K_sideNum','q','T_sidePointNums','len','tk_sideLocs','sk','dimension','b00_startP', 'q_children'});