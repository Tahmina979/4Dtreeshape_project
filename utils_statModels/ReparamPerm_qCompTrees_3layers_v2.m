% Algorithm 2 in paper
function [G,qCompTree_3Ls1p, qCompTree_3Ls2] = ReparamPerm_qCompTrees_3layers_v2(qCompTree_3Ls1, qCompTree_3Ls2, lam_m,lam_s,lam_p)


% Pad minor trees
[qCompTree_3Ls1, qCompTree_3Ls2] = CompatMax_3layers(qCompTree_3Ls1, qCompTree_3Ls2);

%%%%% MATCH SIDES %%%%% 
K1 = qCompTree_3Ls1.K_sideNum;
K2 = qCompTree_3Ls2.K_sideNum;

s1k = qCompTree_3Ls1.sk;
s2k = qCompTree_3Ls2.sk;

% build matching energy and side-branch reparameterization matrices
%% --- Original Method ---

gam_side_all = cell(K1,K2);

E = zeros(K1,K2);

for i=1:K1
    for j=1:K2
        % --- the process here just makes code safe. ---
%         if isempty(qCompTree_3Ls1.q{i})
%             qCompTree_3Ls1.q{i} =[0,  0;
%                                   0,  0;
%                                   0,  0];
%         end
%         
%         if isempty(qCompTree_3Ls2.q{j})
%             qCompTree_3Ls2.q{j}= [0,  0;
%                                   0,  0;
%                                   0,  0];
%         end
        
%         [gam_side_all{i,j}, Eside0] = DPQ_difflen(qCompTree_3Ls1.q{i}, qCompTree_3Ls2.q{j});
        
        % --- compute the overall shape difference between i-th subtree in tree 1 and j-th subtrees of tree 2
        [G, ~] = ReparamPerm_qST_v2(qCompTree_3Ls1.q_children{i}, qCompTree_3Ls2.q_children{j}, lam_m,lam_s,lam_p);
        Eside = sqrt(G.E);
        
        E(i,j) = lam_s*Eside + lam_p*(s1k(i)-s2k(j)).^2;
    end
   
end

% t2 = toc(BB);

% optimize assignment
% a = tic;
[Mvec,Efin] = munkres(E);

% get list of matched pairs and unmatched from munkres output
matched = zeros(0,2);
gam_side = cell(1,0);

for i=1:K1
    matched = [matched; i,Mvec(i)];
end


%% --- Improved Method ---

% gam_side_all = cell(K1,K2);
% E1 = zeros(K1,K2);
% for i=1:K1
%     %fprintf('%d...',i);
%     for j=1:K2
%         [gam_side_all{i,j}, Eside0] = DPQ_difflen(qCompTree_3Ls1.q{i},qCompTree_3Ls2.q{j});
%         [G, ~] = ReparamPerm_qST_new(qCompTree_3Ls1.q_children{i}, qCompTree_3Ls2.q_children{j}, lam_m,lam_s,lam_p);
%         Eside = sqrt(G.E);
%         E1(i,j) = lam_s*Eside + lam_p*(s1k(i)-s2k(j)).^2;
%     end
%     
% end
% 
% % t1 = toc(AA);
% 
% [Mvec1,Efin] = munkres(E1);
% 
% matched = zeros(0,2);
% unmatched1 = [];
% unmatched2 = [];
% gam_side = cell(1,0);
% 
% for i=1:K1
%     if Mvec1(i) ~=0
%         matched = [matched; i, Mvec1(i)];
%         gam_side = [gam_side, gam_side_all{i,Mvec1(i)}];
%     else
%         unmatched1 = [unmatched1, i];
%     end
% end
% 
% for i=1:K2
%     if any(Mvec1 == i) == 0
%         unmatched2 = [unmatched2, i];
%     end
% end

%%


[gam0, Emain] = DPQ_difflen(qCompTree_3Ls1.q0, qCompTree_3Ls2.q0);  % --- DPQ_difflen computes the difference between two q-space branch.

Etotal = lam_m*Emain + Efin;

G = struct('E',Etotal, 'gam0',gam0, 'matched',matched);

%%%%% TRANSFORM qCompTree_3Ls1 %%%%%
%%% q0, t

qCompTree_3Ls1p.q0 = GammaActionQ(qCompTree_3Ls1.q0, G.gam0);
qCompTree_3Ls1p.t_paras = qCompTree_3Ls2.t_paras;

%% s, len0
qCompTree_3Ls1p.s = cumtrapz(qCompTree_3Ls1p.t_paras, sum(qCompTree_3Ls1p.q0.^2,1) );
qCompTree_3Ls1p.len0 = qCompTree_3Ls1p.s(end);

if qCompTree_3Ls1p.len0 == 0
    qCompTree_3Ls1p.s = qCompTree_3Ls1p.t_paras*0.00001;
else
    qCompTree_3Ls1p.s = qCompTree_3Ls1p.s/qCompTree_3Ls1p.len0;
end



m = size(matched,1);
K = m;

%%% qi -- order is [((order of q2)), unmatched1] and q_children
qCompTree_3Ls1p.q = cell(1,K);
qCompTree_3Ls1p.q_children = cell(1, K);
for i=1:m
   
    qCompTree_3Ls1p.q_children{matched(i,2)} = qCompTree_3Ls1.q_children{matched(i,1)}; 
%     qCompTree_3Ls1p.q_children{matched(i,2)}.q0 = qCompTree_3Ls1.q{matched(i,2)};
    qCompTree_3Ls1p.q_children{matched(i,2)}.q0 = qCompTree_3Ls1.q{matched(i,1)};
    qCompTree_3Ls1p.q{matched(i,2)} = qCompTree_3Ls1.q{matched(i,1)};
end


%%% tk, sk
qCompTree_3Ls1p.tk_sideLocs = zeros(1,K);
qCompTree_3Ls1p.sk = zeros(1,K);

% matched,unmatched1 just get updated by gam0
tk1p = interp1(G.gam0, qCompTree_3Ls1p.t_paras, qCompTree_3Ls1.tk_sideLocs);
sk1p = interp1(qCompTree_3Ls1p.t_paras, qCompTree_3Ls1p.s, tk1p);

for i=1:m
    qCompTree_3Ls1p.tk_sideLocs(matched(i,2)) = tk1p(matched(i,1)');
    qCompTree_3Ls1p.sk(matched(i,2)) = sk1p(matched(i,1)');
end

%%% K
qCompTree_3Ls1p.K_sideNum = K;

%%% T0
qCompTree_3Ls1p.T0_pointNum = qCompTree_3Ls2.T0_pointNum;

%%% T
qCompTree_3Ls1p.T_sidePointNums = [];
for i=1:K
    qCompTree_3Ls1p.T_sidePointNums(i) = size(qCompTree_3Ls1p.q{i},2);
end

%%% len
% qCompTree_3Ls1p.len = zeros(1,K);
qCompTree_3Ls1p.len = [];
for k=1:K
    qCompTree_3Ls1p.len(k) = trapz( sum(qCompTree_3Ls1p.q{k}.^2,1) / (qCompTree_3Ls1p.T_sidePointNums(k)-1) );
end

%%% d
qCompTree_3Ls1p.dimension = qCompTree_3Ls1.dimension;

%%% b00
qCompTree_3Ls1p.b00_startP = qCompTree_3Ls1.b00_startP;

% if isfield(qCompTree_3Ls1p, 'q_children') == 0
%     qCompTree_3Ls1p.q_children = cell(1, 0);
% end
qCompTree_3Ls1p = orderfields(qCompTree_3Ls1p, ...
    {'t_paras','s','q0','T0_pointNum','len0','K_sideNum','q','T_sidePointNums','len','tk_sideLocs','sk','dimension','b00_startP','q_children'});

%% --- Do the matching of layer 2 and layer 3 ---

for i = 1: numel(qCompTree_3Ls1p.q_children)
    
    % --- 
    [G1, qCompTree_3Ls1p.q_children{i}, qCompTree_3Ls2.q_children{i}] = ReparamPerm_qST_v2(qCompTree_3Ls1p.q_children{i}, qCompTree_3Ls2.q_children{i},lam_m,lam_s,lam_p);
        
%     if numel(qCompTree_3Ls1p.q_children_p{i}.q) > numel(qCompTree_3Ls2p.q_children{i}.q)
%         qCompTree_3Ls2p.q_children_p{i} = wg_root_AddZeroBranches(qCompTree_3Ls1p.q_children_p{i}, qCompTree_3Ls2p.q_children{i});
%     else
%         qCompTree_3Ls2p.q_children_p{i} = qCompTree_3Ls2p.q_children{i};
%     end
end


% --- here we should make the layer 2 branches's number equal for the two trees ---

% if isfield(qCompTree_3Ls2, 'q_children') == 0
%     qCompTree_3Ls2.q_children = cell(1, 0);
% end
% 
% 
% if numel(qCompTree_3Ls1p.q_children) > numel(qCompTree_3Ls2.q_children)
%     qCompTree_3Ls2p = wg_root_AddZeroBranches(qCompTree_3Ls1p, qCompTree_3Ls2);
% else
%     qCompTree_3Ls2p = qCompTree_3Ls2;
% end
% 
% 
% for i = 1: numel(qCompTree_3Ls1p.q_children)
%     
%     % --- 
%     [G1, qCompTree_3Ls1p.q_children_p{i}, qCompTree_3Ls2p.q_children{i}] = ReparamPerm_qST_wg(qCompTree_3Ls1p.q_children{i}, qCompTree_3Ls2p.q_children{i},lam_m,lam_s,lam_p);
%         
%     if numel(qCompTree_3Ls1p.q_children_p{i}.q) > numel(qCompTree_3Ls2p.q_children{i}.q)
%         qCompTree_3Ls2p.q_children_p{i} = wg_root_AddZeroBranches(qCompTree_3Ls1p.q_children_p{i}, qCompTree_3Ls2p.q_children{i});
%     else
%         qCompTree_3Ls2p.q_children_p{i} = qCompTree_3Ls2p.q_children{i};
%     end
% end
% 
% if isfield(qCompTree_3Ls1p, 'q_children_p') == 1
%     qCompTree_3Ls1p.q_children = qCompTree_3Ls1p.q_children_p;
% end
% 
% if isfield(qCompTree_3Ls2p, 'q_children_p') == 1
%     qCompTree_3Ls2p.q_children = qCompTree_3Ls2p.q_children_p;
% end




end

