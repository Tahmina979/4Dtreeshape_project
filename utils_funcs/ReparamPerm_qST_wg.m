% Algorithm 2 in paper
function [G,qST1p, qST2] = ReparamPerm_qST_wg(qST1,qST2, lam_m,lam_s,lam_p)
% function [matched,Etotal,E] = ReparamPerm_qST(qST1,qST2, lam_m,lam_s,lam_p)

%%%%% MATCH SIDES %%%%% 
K1 = qST1.K_sideNum;
K2 = qST2.K_sideNum;

s1k = qST1.sk;
s2k = qST2.sk;

% build matching energy and side-branch reparameterization matrices
% --- Case1: The matrix E1 is (K1+ K2)^2 ---
A1 = tic;
A2 = tic;
gam_side_all = cell(K1,K2);
% E1 = zeros(K1, K2);
E = zeros(K1+K2,K1+K2);
for i=1:K1
    %fprintf('%d...',i);
    for j=1:K2
        % --- the process here just makes code safe. ---
        if isempty(qST1.q{i})
            qST1.q{i} = [0,0;0,0;0,0];
        end
        if isempty(qST2.q{j})
            qST2.q{j}= [0,0;0,0;0,0];
        end
        
        [gam_side_all{i,j}, Eside] = DPQ_difflen(qST1.q{i},qST2.q{j});
        E(i,j) = lam_s*Eside + lam_p*(s1k(i)-s2k(j)).^2;
    end
    E(i,K2+1:end) = lam_s*qST1.len(i);
end

for j=1:K2
    E(K1+1:end,j) = lam_s*qST2.len(j);
end

t1_2 = toc(A2);     % --- t1_2 is the just the time of computing DM ---

% optimize assignment

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


%fprintf('\nAlign main...\n');

[gam0, Emain] = DPQ_difflen(qST1.q0,qST2.q0);  % --- DPQ_difflen computes the difference between two q-space branch.
Etotal = lam_m*Emain + Efin;

G = struct('E',Etotal, 'gam0',gam0, 'gam',{gam_side}, ...
    'matched',matched, 'unmatched1',unmatched1, 'unmatched2',unmatched2);

%%%%% TRANSFORM qST1 %%%%%
%%% q0, t
qST1p.q0 = GammaActionQ(qST1.q0, G.gam0);
qST1p.t_paras = qST2.t_paras;

%%% s, len0
if norm(sum(qST1p.q0.^2,1)) == 0
    qST1p.s = qST1p.t_paras *0;
else
    qST1p.s = cumtrapz(qST1p.t_paras, sum(qST1p.q0.^2,1) );
end

qST1p.len0 = qST1p.s(end);

% qST1p.s = qST1p.s/(qST1p.len0+0.000001)

if qST1p.len0 ~= 0
    qST1p.s = qST1p.s/qST1p.len0;
else
    qST1p.s = rand(size(qST1p.s, 1), size(qST1p.s, 2))*0.0000001;
end


m = size(matched,1);
um1= numel(unmatched1);
um2= numel(unmatched2);
K = m + um2 + um1;

%%% qi -- order is [((order of q2)), unmatched1]
qST1p.q = cell(1,K);
for i=1:m
    qST1p.q{matched(i,2)} = GammaActionQ(qST1.q{matched(i,1)}, G.gam{i});
    
end

for i=1:um2
    qST1p.q{unmatched2(i)} = zeros(qST1.dimension,qST2.T_sidePointNums(unmatched2(i)));
end

for i=1:um1
    qST1p.q{m+um2+i} = qST1.q{unmatched1(i)};
end

%%% tk, sk
qST1p.tk_sideLocs = zeros(1,K);
qST1p.sk = zeros(1,K);

% matched,unmatched1 just get updated by gam0
if isempty(qST1.tk_sideLocs)
    tk1p = qST1.tk_sideLocs;
else
    tk1p = interp1(G.gam0, qST1p.t_paras, qST1.tk_sideLocs);
end

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
% if isempty(unmatched2) 
%     qST1p.tk(unmatched2) = [];
% else
%     qST1p.tk(unmatched2) = interp1(qST1p.s, qST1p.t, qST1p.sk(unmatched2));
% end

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
qST1p.len = zeros(1,K);
for k=1:K
    qST1p.len(k) = trapz( sum(qST1p.q{k}.^2,1) / (qST1p.T_sidePointNums(k)-1) );
end

%%% d
qST1p.dimension = qST1.dimension;

%%% b00
qST1p.b00_startP = qST1.b00_startP;

% --- post processing of qST2 to make it safe ---
for i=1: qST2.K_sideNum
    qST2.T_sidePointNum(i) = size(qST2.q{i},2);
end


if isfield(qST1p, 'q_children') == 0
    qST1p = orderfields(qST1p, ...
        {'t_paras','s','q0','T0_pointNum','len0','K_sideNum','q','T_sidePointNums','len','tk_sideLocs','sk','dimension','b00_startP'});
else
    qST1p = orderfields(qST1p, ...
        {'t_paras','s','q0','T0_pointNum','len0','K_sideNum','q','T_sidePointNums','len','tk_sideLocs','sk','dimension','b00_startP','q_children'});
end

end

