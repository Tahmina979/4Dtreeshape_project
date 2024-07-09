function q_zero = makeOneVirtualTree_3layers(q_children)


q_zero.t_paras = q_children.t_paras;
q_zero.s = q_children.s;
q_zero.q0 = zeros(size(q_children.q0, 1), size(q_children.q0, 2));
q_zero.T0_pointNum = q_children.T0_pointNum;
q_zero.len0 = 0;
% q_zero.len0 = q_children.len0;
q_zero.K_sideNum = q_children.K_sideNum;
q_zero.q = q_children.q;


for i = 1: numel(q_zero.q)
    q_zero.q{i} = zeros(size(q_children.q{i}, 1), size(q_children.q{i}, 2));
    q_zero.q_children{i} = makeZeroST(q_children.q_children{i});
end

% q_zero.T = zeros(size(q_children.T, 1), size(q_children.T, 2));
q_zero.T_sidePointNums = q_children.T_sidePointNums;
q_zero.len = zeros(size(q_children.len, 1), size(q_children.len, 2));
% q_zero.tk = zeros(size(q_children.tk, 1), size(q_children.tk, 2));
% q_zero.sk = zeros(size(q_children.sk, 1), size(q_children.sk, 2));
q_zero.tk_sideLocs = q_children.tk_sideLocs;
q_zero.sk = q_children.sk;
q_zero.dimension = q_children.dimension;
q_zero.b00_startP = zeros(size(q_children.b00_startP, 1), size(q_children.b00_startP, 2));

% q_zero.func_t_parasAndRad = [];



end