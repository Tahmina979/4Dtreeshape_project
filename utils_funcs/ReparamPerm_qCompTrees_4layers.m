% Algorithm 2 in paper
function [G,qCompTree_4Ls1_p, qCompTree_4Ls2_p] = ReparamPerm_qCompTrees_4layers(qCompTree_4Ls1,qCompTree_4Ls2, lam_m,lam_s,lam_p)

% Pad minor trees
% [qCompTree_4Ls1, qCompTree_4Ls2] = CompatMax_4layers(qCompTree_4Ls1, qCompTree_4Ls2);


%%%%% MATCH SIDES %%%%% 
K1 = qCompTree_4Ls1.K_sideNum;
K2 = qCompTree_4Ls2.K_sideNum;

s1k = qCompTree_4Ls1.sk;
s2k = qCompTree_4Ls2.sk;

% build matching energy and side-branch reparameterization matrices
%% Ducun's method

% --- Energy matrix computation ---
gam_side_all = cell(K1,K2);

E = zeros(K1+K2,K1+K2);

for i=1:K1
    %fprintf('%d...',i);
    for j=1:K2
        % --- the process here just makes code safe. ---
        if isempty(qCompTree_4Ls1.q{i})
            qCompTree_4Ls1.q{i} = [0,0;0,0;0,0];
        end
        
        if isempty(qCompTree_4Ls2.q{j})
            qCompTree_4Ls2.q{j}= [0,0;0,0;0,0];
        end
        
%         [gam_side_all{i,j}, Eside0] = DPQ_difflen(qCompTree_4Ls1.q{i}, qCompTree_4Ls2.q{j});
        [G, ~] = ReparamPerm_qCompTrees_3layers(qCompTree_4Ls1.q_children{i}, ... 
                                                qCompTree_4Ls2.q_children{j}, lam_m,lam_s,lam_p);
        Eside = sqrt(G.E);
        
        E(i,j) = lam_s*Eside + lam_p*(s1k(i)-s2k(j)).^2;
    end
   
%     [G1] = ReparamPerm_qST_complex_forZero(qCompTree_4Ls1.q_children{i}, q_zero, lam_m,lam_s,lam_p);
    
    % --- 
    Sum_len = qCompTree_4Ls1.q_children{i}.len0;
    
    for t = 1: numel(qCompTree_4Ls1.q_children{i}.q_children)
        Sum_len = Sum_len + qCompTree_4Ls1.q_children{i}.q_children{t}.len0 + sum(qCompTree_4Ls1.q_children{i}.q_children{t}.len);
    end
    
    
    E(i,K2+1:end) = lam_s* Sum_len;
end

for j=1:K2
%     clear q_zero
%     q_zero = makeZeroComplexTree(qCompTree_4Ls2.q_children{j});
%     [G2] = ReparamPerm_qST_complex_forZero(q_zero, qCompTree_4Ls2.q_children{j}, lam_m,lam_s,lam_p);
%     
    subTree = qCompTree_4Ls2.q_children{j};
    Sum_len = subTree.len0;
    
    for t = 1: numel(subTree.q_children)
        Sum_len = Sum_len + subTree.q_children{t}.len0 + sum(subTree.q_children{t}.len);
    end
    E(K1+1:end,j) = lam_s* Sum_len;
end
% t2 = toc(BB);

%% optimize assignment, Compute correspondence
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

%fprintf('\nAlign main...\n');

[gam0, Emain] = DPQ_difflen(qCompTree_4Ls1.q0, qCompTree_4Ls2.q0);  % --- DPQ_difflen computes the difference between two q-space branch.

Etotal = lam_m*Emain + Efin;

G = struct( 'E',Etotal, ...
            'gam0',gam0, ...
            'gam',{gam_side}, ...
            'matched',matched, ...
            'unmatched1',unmatched1,...
            'unmatched2',unmatched2);
        
%%
%%%%% TRANSFORM qCompTree_4Ls1 %%%%%
%%% q0, t
qCompTree_4Ls1_p.q0 = GammaActionQ(qCompTree_4Ls1.q0, G.gam0);
qCompTree_4Ls1_p.t_paras = qCompTree_4Ls2.t_paras;

%% s, len0
qCompTree_4Ls1_p.s = cumtrapz(qCompTree_4Ls1_p.t_paras, sum(qCompTree_4Ls1_p.q0.^2,1) );
qCompTree_4Ls1_p.len0 = qCompTree_4Ls1_p.s(end);

if qCompTree_4Ls1_p.len0 == 0
    qCompTree_4Ls1_p.s = qCompTree_4Ls1_p.t_paras*0.00001;
else
    qCompTree_4Ls1_p.s = qCompTree_4Ls1_p.s/qCompTree_4Ls1_p.len0;
end



m = size(matched,1);
um1= numel(unmatched1);
um2= numel(unmatched2);
K = m + um2 + um1;

%%% qi -- order is [((order of q2)), unmatched1]
qCompTree_4Ls1_p.q = cell(1,K);
qCompTree_4Ls2_p.q = cell(1,K);

for i=1:m
    
    qCompTree_4Ls1_p.q{matched(i,2)} = GammaActionQ(qCompTree_4Ls1.q{matched(i,1)}, G.gam{i});
    qCompTree_4Ls1_p.q_children{matched(i,2)} = qCompTree_4Ls1.q_children{matched(i,1)};
    
    qCompTree_4Ls1_p.q_children{matched(i,2)}.q0 = qCompTree_4Ls1_p.q{matched(i,2)};
%     qCompTree_4Ls2_p.q{matched(i,2)} = qCompTree_4Ls2.q{matched(i,2)};
%     qCompTree_4Ls1_p.q_children{matched(i,2)} = qCompTree_4Ls1.q_children{matched(i,2)};
end

for i=1:um2
    
    qCompTree_4Ls1_p.q{unmatched2(i)} = zeros(qCompTree_4Ls1.dimension, qCompTree_4Ls2.T_sidePointNums(unmatched2(i)));
    qCompTree_4Ls1_p.q_children{unmatched2(i)} = makeZeroComplexTree(qCompTree_4Ls2.q_children{unmatched2(i)});
    qCompTree_4Ls1_p.q_children{unmatched2(i)}.q0 = qCompTree_4Ls1_p.q{unmatched2(i)};
end

for i=1:um1
    
    qCompTree_4Ls1_p.q{m+um2+i} = qCompTree_4Ls1.q{unmatched1(i)};
    qCompTree_4Ls1_p.q_children{m+um2+i} = qCompTree_4Ls1.q_children{unmatched1(i)};
    qCompTree_4Ls1_p.q_children{m+um2+i}.q0 = qCompTree_4Ls1_p.q{m+um2+i};
%     qCompTree_4Ls2_p.q{m+um2+i} = zeros(qCompTree_4Ls1.d, qCompTree_4Ls1.T(unmatched1(i)));
%     qCompTree_4Ls2_p.q_children{m+um2+i} = wg_root_makeZeroST(qCompTree_4Ls1.q_children{unmatched1(i)});
%     
end

%%% tk, sk
qCompTree_4Ls1_p.tk_sideLocs = zeros(1,K);
qCompTree_4Ls1_p.sk = zeros(1,K);

% matched,unmatched1 just get updated by gam0
tk1p = interp1(G.gam0, qCompTree_4Ls1_p.t_paras, qCompTree_4Ls1.tk_sideLocs);
sk1p = interp1(qCompTree_4Ls1_p.t_paras, qCompTree_4Ls1_p.s, tk1p);

for i=1:m
    qCompTree_4Ls1_p.tk_sideLocs(matched(i,2)) = tk1p(matched(i,1)');
    qCompTree_4Ls1_p.sk(matched(i,2)) = sk1p(matched(i,1)');
end

for i=1:um1
    qCompTree_4Ls1_p.tk_sideLocs(m+um2+i) = tk1p(unmatched1(i));
    qCompTree_4Ls1_p.sk(m+um2+i) = sk1p(unmatched1(i));
end

% unmatched2 get arclength position from tree 2
qCompTree_4Ls1_p.sk(unmatched2) = s2k(unmatched2);
qCompTree_4Ls1_p.tk_sideLocs(unmatched2) = interp1(qCompTree_4Ls1_p.s, qCompTree_4Ls1_p.t_paras, qCompTree_4Ls1_p.sk(unmatched2));

%%% K
qCompTree_4Ls1_p.K_sideNum = K;

%%% T0
qCompTree_4Ls1_p.T0_pointNum = qCompTree_4Ls2.T0_pointNum;

%%% T
qCompTree_4Ls1_p.T_sidePointNums = [];
for i=1:K
    qCompTree_4Ls1_p.T_sidePointNums(i) = size(qCompTree_4Ls1_p.q{i},2);
end

%%% len
% qCompTree_4Ls1_p.len = zeros(1,K);
qCompTree_4Ls1_p.len = [];
for k=1:K
    qCompTree_4Ls1_p.len(k) = trapz( sum(qCompTree_4Ls1_p.q{k}.^2,1) / (qCompTree_4Ls1_p.T_sidePointNums(k)-1) );
end

%%% d
qCompTree_4Ls1_p.dimension = qCompTree_4Ls1.dimension;

%%% b00
qCompTree_4Ls1_p.b00_startP = qCompTree_4Ls1.b00_startP;


qCompTree_4Ls1_p = orderfields(qCompTree_4Ls1_p, ...
    {'t_paras','s','q0','T0_pointNum','len0','K_sideNum','q','T_sidePointNums','len','tk_sideLocs','sk','dimension','b00_startP','q_children'});

%% --- Do the matching of layer 2 and layer 3 ---

% --- here we should make the layer 2 branches equal for the two trees ---
if numel(qCompTree_4Ls1_p.q_children) > numel(qCompTree_4Ls2.q_children)
%     qCompTree_4Ls2_p = wg_root_AddZeroBranches(qCompTree_4Ls1_p, qCompTree_4Ls2);
    qCompTree_4Ls2_p = AddVitualTrees_3layers(qCompTree_4Ls1_p, qCompTree_4Ls2);

else
    qCompTree_4Ls2_p = qCompTree_4Ls2;
end


% ---
for i = 1: numel(qCompTree_4Ls1_p.q_children)
    
    % --- here it is better that numel(qCompTree_4Ls1_p.q_children{i}.q) is less than
    % numel(qCompTree_4Ls2_p.q_children{i}.q) .
    [G1, qCompTree_4Ls1_p.q_children_p{i}, qCompTree_4Ls2_p.q_children_p{i}] = ...
                            ReparamPerm_qCompTrees_3layers(qCompTree_4Ls1_p.q_children{i}, qCompTree_4Ls2_p.q_children{i}, lam_m,lam_s,lam_p);
                            
%     [G, qCompTree_4Ls1_p.q_children_p{i}] = ReparamPerm_qST(qCompTree_4Ls1_p.q_children{i}, qCompTree_4Ls2_p.q_children{j}, lam_m,lam_s,lam_p);
        
%     if numel(qCompTree_4Ls1_p.q_children_p{i}.q) > numel(qCompTree_4Ls2_p.q_children{i}.q)
%         qCompTree_4Ls2_p.q_children_p{i} = AddVitualTrees_3layers(qCompTree_4Ls1_p.q_children_p{i}, qCompTree_4Ls2_p.q_children{i});
%     else
%         qCompTree_4Ls2_p.q_children_p{i} = qCompTree_4Ls2_p.q_children{i};
%     end
end

qCompTree_4Ls1_p.q_children = qCompTree_4Ls1_p.q_children_p;
qCompTree_4Ls2_p.q_children = qCompTree_4Ls2_p.q_children_p;



end

