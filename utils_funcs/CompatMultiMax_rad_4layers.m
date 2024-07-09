function [Qs] = CompatMultiMax_rad_4layers(Qs)


%%% main branch
T0_pointNum_all = cellfun(@(A)A.T0_pointNum, Qs);
T0_pointNum = max(T0_pointNum_all);
t = linspace(0,1,T0_pointNum);
q0s = [];
beta0_rads =[];
% q0 and beta0_rad

[q0s, beta0_rads] = cellfun(@(Q) ReCalc_q0(Q, T0_pointNum,t), Qs, 'UniformOutput', true);

%%% side branches
% Km = min(Q1.K_sideNum,Q2.K_sideNum);
% K = max(Q1.K_sideNum,Q2.K_sideNum);

K_all = cellfun(@(A)A.K_sideNum, Qs);
sk_all = cellfun(@(A)A.sk, Qs, 'UniformOutput', false);
beta_rad_all = cellfun(@(A)A.beta_rad, Qs, 'UniformOutput', false);
q_children_3Ls_all = cellfun(@(A)A.q_children, Qs, 'UniformOutput', false);
[K_max, max_id] = max(K_all);

% sk
sk_max = sk_all{max_id};
beta_rad_max = beta_rad_all{max_id};
q_children_3Ls_max = q_children_3Ls_all{max_id};

[sk_s, beta_rads] = cellfun(@(A) ReCalc_sk_beta_rad(A, K_max, sk_max, beta_rad_max), ...
                            Qs, 'UniformOutput', false);

% 3layers sub trees
q_children_3Ls = cellfun(@(A) ReCalc_3LsSubtrees(A, K_max, q_children_3Ls_max),...
                        Qs, 'UniformOutput', false);




% make trees
% Qs = cell(1, length(Qs));
for i=1: length(Qs)
    Qs{i} = make_qCompTree_rad_4layers(q0s{i}, q_children_3Ls{i}, sk_s{i}, Qs{i}.b00_startP, ...
                                              beta0_rads{i}, beta_rads{i});
end


%% 3Layers Subtrees
treeNum= length(Qs);

for i=1: length(Qs{1}.q_children)
    q_children_all = cellfun(@(A)A.q_children{i}, Qs, 'UniformOutput', false);
    [q_children_all] = CompatMultiMax_rad_3layers(q_children_all);
    
    
    for t=1: treeNum
        Qs{t}.q_children{i} = q_children_all{t};
        [Qs{t},~] = ResetqCompTree_q_and_beta_rad(Qs{t},Qs{t});
    end
    
end



% 
% % 3layers sub trees
% q_children_3Ls1 = cell(1,K);
% q_children_3Ls2 = cell(1,K);
% 
% for k=1:Km
%     q_children_3Ls1{k} = Q1.q_children{k};
%     q_children_3Ls2{k} = Q2.q_children{k};
% end
% 
% if Q1.K_sideNum<K
%     for k=Km+1:K
%         q_children_3Ls2{k} = Q2.q_children{k};
%         q_children_3Ls1{k} = makeOneVirtualTree_rad_3layers(q_children_3Ls2{k});
%     end
% elseif Q2.K_sideNum<K
%     for k=Km+1:K
%         q_children_3Ls1{k} = Q1.q_children{k};
%         q_children_3Ls2{k} = makeOneVirtualTree_rad_3layers(q_children_3Ls1{k});
%     end
% end
% 
% % make return values
% Q1p = make_qCompTree_rad_4layers(q0_1, q_children_3Ls1, sk_1, Q1.b00_startP, ...
%                                 beta0_rad_1, beta_rad_1);
%                             
% Q2p = make_qCompTree_rad_4layers(q0_2, q_children_3Ls2, sk_2, Q2.b00_startP, ...
%                                 beta0_rad_2, beta_rad_2);

end


function [q0, beta0_rad] =ReCalc_q0(Q, T0_pointNum, t)

if Q.T0_pointNum<T0_pointNum
    q0 = interp1(Q.t_paras,Q.q0', t)'; 
    beta0_rad = interp1(Q.t_paras, Q.beta0_rad, t);
else
    q0 = Q.q0;
    beta0_rad = Q.beta0_rad;
end

q0 = {q0};
beta0_rad = {beta0_rad};
end

function [sk, beta_rad]= ReCalc_sk_beta_rad(Q, K_max, sk_max, beta_rad_max)

sk = Q.sk;
beta_rad = Q.beta_rad;
K = Q.K_sideNum;
if K < K_max
    sk = [sk, sk_max(K+1:K_max)];
    beta_rad = [beta_rad, beta_rad_max(K+1:K_max)];
end

end

function [q_children_3Ls]= ReCalc_3LsSubtrees(Q, K_max, q_children_3Ls_max)

K = Q.K_sideNum;
q_children_3Ls = cell(1,K);
for i=1:K
    q_children_3Ls{i} = Q.q_children{i};
end

if K < K_max
    for i=K+1:K_max
        q_children_3Ls{i} = makeOneVirtualTree_rad_3layers(q_children_3Ls_max{i});
    end
end


end




