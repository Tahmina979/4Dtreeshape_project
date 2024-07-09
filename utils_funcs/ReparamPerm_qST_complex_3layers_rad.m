% Algorithm 2 in paper
function [G,qST1p, qST2p] = ReparamPerm_qST_complex_3layers_rad(qST1,qST2, lam_m,lam_s,lam_p)

%%%%% MATCH SIDES %%%%% 
K1 = qST1.K_sideNum;
K2 = qST2.K_sideNum;

s1k = qST1.sk;
s2k = qST2.sk;

% build matching energy and side-branch reparameterization matrices
%% --- Original Method ---

gam_side_all = cell(K1,K2);

E = zeros(K1+K2,K1+K2);

for i=1:K1
     
    for j=1:K2
        % --- the process here just makes code safe. ---
        if isempty(qST1.q{i})
            qST1.q{i} = [0,0;0,0;0,0; 0.1,0.1];
        end
        
        if isempty(qST2.q{j})
            qST2.q{j}= [0,0;0,0;0,0;0.1,0.1];
        end
        
        [gam_side_all{i,j}, Eside0] = DPQ_difflen(qST1.q{i}, qST2.q{j});
        
        % --- compute the overall shape difference between i-th subtree in tree 1 and j-th subtrees of tree 2
        [G, ~] = ReparamPerm_qST_rad(qST1.q_children{i}, qST2.q_children{j}, lam_m,lam_s,lam_p);
        Eside = sqrt(G.E);
        
        E(i,j) = lam_s*Eside + lam_p*(s1k(i)-s2k(j)).^2;
    end
    
    q_zero = makeZeroST(qST1.q_children{i});
    [G1] = ReparamPerm_qST_forZero_rad(qST1.q_children{i}, q_zero, lam_m,lam_s,lam_p);
    E(i,K2+1:end) = lam_s* sqrt(G1.E);
end

for j=1:K2
    clear q_zero
    q_zero = makeZeroST(qST2.q_children{j});
    [G2] = ReparamPerm_qST_forZero_rad(q_zero, qST2.q_children{j}, lam_m,lam_s,lam_p);
    
    E(K1+1:end,j) = lam_s* sqrt(G2.E);
end
% t2 = toc(BB);

% optimize assignment
% a = tic;
[Mvec,Efin] = munkres(E);

% get list of matched pairs and unmatched from munkres output
matched = zeros(0,2);
unmatched1 = [];
unmatched2 = [];
gam_side = cell(1,0);

for i=1:K1
    if Mvec(i)<=K2
        matched = [matched; i,Mvec(i)];
        gam_side = [gam_side, gam_side_all{i,Mvec(i)}];
    else
        unmatched1 = [unmatched1, i];
    end
end

for i= K1+(1:K2)
    if Mvec(i)<=K2
        unmatched2 = [unmatched2, Mvec(i)];
    end
end


%%


[gam0, Emain] = DPQ_difflen(qST1.q0, qST2.q0);  % --- DPQ_difflen computes the difference between two q-space branch.

Etotal = lam_m*Emain + Efin;

G = struct('E',Etotal, 'gam0',gam0, 'gam',{gam_side}, ...
    'matched',matched, 'unmatched1',unmatched1, 'unmatched2',unmatched2);

%%%%% TRANSFORM qST1 %%%%%
%%% q0, t

qST1p.q0 = GammaActionQ(qST1.q0, G.gam0);
qST1p.t_paras = qST2.t_paras;

%% s, len0
qST1p.s = cumtrapz(qST1p.t_paras, sum(qST1p.q0.^2,1) );
qST1p.len0 = qST1p.s(end);

if qST1p.len0 == 0
    qST1p.s = qST1p.t_paras*0.00001;
else
    qST1p.s = qST1p.s/qST1p.len0;
end

% --- Extra added it here ---
qST1p.s = qST1p.t_paras;

m = size(matched,1);
um1= numel(unmatched1);
um2= numel(unmatched2);
K = m + um2 + um1;

%%% qi -- order is [((order of q2)), unmatched1]
qST1p.q = cell(1,K);
qST2p.q = cell(1,K);

for i=1:m
    
    qST1p.q{matched(i,2)} = GammaActionQ(qST1.q{matched(i,1)}, G.gam{i});
    qST1p.q_children{matched(i,2)} = qST1.q_children{matched(i,1)};
    
    qST1p.q_children{matched(i,2)}.q0 = qST1p.q{matched(i,2)};
%     qST2p.q{matched(i,2)} = qST2.q{matched(i,2)};
%     qST1p.q_children{matched(i,2)} = qST1.q_children{matched(i,2)};
end

for i=1:um2
    
    qST1p.q{unmatched2(i)} = zeros(qST1.dimension, qST2.T_sidePointNums(unmatched2(i)));
    qST1p.q_children{unmatched2(i)} = makeZeroST(qST2.q_children{unmatched2(i)});
    qST1p.q_children{unmatched2(i)}.q0 = qST1p.q{unmatched2(i)};
end

for i=1:um1
    
    qST1p.q{m+um2+i} = qST1.q{unmatched1(i)};
    qST1p.q_children{m+um2+i} = qST1.q_children{unmatched1(i)};
    qST1p.q_children{m+um2+i}.q0 = qST1p.q{m+um2+i};
    
end

%%% tk, sk
qST1p.tk_sideLocs = zeros(1,K);
qST1p.sk = zeros(1,K);

% matched,unmatched1 just get updated by gam0
tk1p = interp1(G.gam0, qST1p.t_paras, qST1.tk_sideLocs);
tk1p = qST1.tk_sideLocs;

sk1p = interp1(qST1p.t_paras, qST1p.s, tk1p);

for i=1:m
    qST1p.tk_sideLocs(matched(i,2)) = tk1p(matched(i,1)');
    qST1p.sk(matched(i,2)) = sk1p(matched(i,1)');
end

for i=1:um1
    qST1p.tk_sideLocs(m+um2+i) = tk1p(unmatched1(i));
    qST1p.sk(m+um2+i) = sk1p(unmatched1(i));
end

% unmatched2 get arclength position from tree 2
qST1p.sk(unmatched2) = s2k(unmatched2);
qST1p.tk_sideLocs(unmatched2) = interp1(qST1p.s, qST1p.t_paras, qST1p.sk(unmatched2));

%%% K
qST1p.K_sideNum = K;

%%% T0
qST1p.T0_pointNum = qST2.T0_pointNum;

%%% T
qST1p.T_sidePointNums = [];
for i=1:K
    qST1p.T_sidePointNums(i) = size(qST1p.q{i},2);
end

%%% len
% qST1p.len = zeros(1,K);
qST1p.len = [];
for k=1:K
    qST1p.len(k) = trapz( sum(qST1p.q{k}.^2,1) / (qST1p.T_sidePointNums(k)-1) );
end

%%% d
qST1p.dimension = qST1.dimension;

%%% b00
qST1p.b00_startP = qST1.b00_startP;

if isfield(qST1p, 'q_children') == 0
    qST1p.q_children = cell(1, 0);
end

qST1p.func_t_parasAndRad = qST1.func_t_parasAndRad;
qST1p.data_t_parasAndRad = qST1.data_t_parasAndRad;

if isempty(qST1p.func_t_parasAndRad) == 0
    qST1p.q0(4, :) = ppval(qST1p.func_t_parasAndRad, qST1p.t_paras);
end

if isempty(qST1p.data_t_parasAndRad) == 0
    qST1p.q0(4, :) = interp1(qST1p.data_t_parasAndRad(1, :), qST1p.data_t_parasAndRad(2, :), qST1p.t_paras);
end


qST1p = orderfields(qST1p, ...
    {'t_paras','s','q0','T0_pointNum','len0','K_sideNum','q','T_sidePointNums','len','tk_sideLocs','sk','dimension','b00_startP','q_children','func_t_parasAndRad', 'data_t_parasAndRad'});

%% --- Do the matching of layer 2 and layer 3 ---

% --- here we should make the layer 2 branches's number equal for the two trees ---

if isfield(qST2, 'q_children') == 0
    qST2.q_children = cell(1, 0);
end


if numel(qST1p.q_children) > numel(qST2.q_children)
    qST2p = wg_root_AddZeroBranches(qST1p, qST2);
else
    qST2p = qST2;
end


% --- re-compute the radius ---
qST2p.func_t_parasAndRad = qST2.func_t_parasAndRad;

% --- Here we do this for safe ---
if isfield(qST2, 'data_t_parasAndRad') == 0
    qST2.data_t_parasAndRad = [];
end

qST2p.data_t_parasAndRad = qST2.data_t_parasAndRad;

if isempty(qST2p.func_t_parasAndRad) == 0
    qST2p.q0(4, :) = ppval(qST2p.func_t_parasAndRad, qST2p.t_paras);
end

if isempty(qST2p.data_t_parasAndRad) == 0
    qST2p.q0(4, :) = interp1(qST2p.data_t_parasAndRad(1, :), qST2p.data_t_parasAndRad(2, :), qST2p.t_paras);
end


for i = 1: numel(qST1p.q_children)
    
    % --- 
    [G1, qST1p.q_children_p{i}, qST2p.q_children{i}] = ReparamPerm_qST_wg_rad(qST1p.q_children{i}, qST2p.q_children{i},lam_m,lam_s,lam_p);
        
    if numel(qST1p.q_children_p{i}.q) > numel(qST2p.q_children{i}.q)
        qST2p.q_children_p{i} = wg_root_AddZeroBranches(qST1p.q_children_p{i}, qST2p.q_children{i});
    else
        qST2p.q_children_p{i} = qST2p.q_children{i};
        qST2p.q_children_p{i}.func_t_parasAndRad = qST2p.q_children{i}.func_t_parasAndRad;
        qST2p.q_children_p{i}.data_t_parasAndRad = qST2p.q_children{i}.data_t_parasAndRad;
        
        
        if isempty(qST2p.q_children_p{i}.func_t_parasAndRad) == 0
            qST2p.q_children_p{i}.q0(4, :) = ppval(qST2p.q_children_p{i}.func_t_parasAndRad, qST2p.q_children_p{i}.t_paras);
        end
        
        if isempty(qST2p.q_children_p{i}.data_t_parasAndRad) == 0
            qST2p.q_children_p{i}.q0(4, :) = interp1(qST2p.q_children_p{i}.data_t_parasAndRad(1, :), qST2p.q_children_p{i}.data_t_parasAndRad(2, :), ......
                                                                       qST2p.q_children_p{i}.t_paras);
        end
        
    end
    
    % --- recompute the radius ---
    
% 
%     if isempty(qST2p.q_children_p{i}.func_t_parasAndRad) == 0
%         qST2p.q_children_p{i}.q0(4, :) = ppval(qST2p.q_children_p{i}.func_t_parasAndRad, qST2p.q_children_p{i}.t_paras);
%     end
    
    
end

if isfield(qST1p, 'q_children_p') == 1
    qST1p.q_children = qST1p.q_children_p;
end

if isfield(qST2p, 'q_children_p') == 1
    qST2p.q_children = qST2p.q_children_p;
end




end

