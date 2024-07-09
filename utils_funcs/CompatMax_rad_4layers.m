function [Q1p,Q2p] = CompatMax_rad_4layers(Q1, Q2)


%%% main branch
T0_pointNum = max(Q1.T0_pointNum, Q2.T0_pointNum);
t = linspace(0,1,T0_pointNum);

% q0
if Q1.T0_pointNum<T0_pointNum
    q0_1 = interp1(Q1.t_paras,Q1.q0', t)'; 
    beta0_rad_1 = interp1(Q1.t_paras, Q1.beta0_rad, t);
    
    q0_2 = Q2.q0;
    beta0_rad_2 = Q2.beta0_rad;
elseif Q2.T0_pointNum<T0_pointNum
    q0_1 = Q1.q0;
    beta0_rad_1 = Q1.beta0_rad;
    
    q0_2 = interp1(Q2.t_paras,Q2.q0', t)';
    beta0_rad_2 = interp1(Q2.t_paras, Q2.beta0_rad, t);
else
    q0_1 = Q1.q0; 
    beta0_rad_1 = Q1.beta0_rad;
    q0_2 = Q2.q0; 
    beta0_rad_2 = Q2.beta0_rad;
end

%%% side branches
Km = min(Q1.K_sideNum,Q2.K_sideNum);
K = max(Q1.K_sideNum,Q2.K_sideNum);

% sk
sk_1 = Q1.sk;
sk_2 = Q2.sk;
if Q1.K_sideNum<K
    sk_1 = [sk_1, sk_2(Km+1:K)];
elseif Q2.K_sideNum<K
    sk_2 = [sk_2, sk_1(Km+1:K)];
end

% beta rad
beta_rad_1 = Q1.beta_rad;
beta_rad_2 = Q2.beta_rad;
if Q1.K_sideNum<K
    beta_rad_1 = [beta_rad_1, beta_rad_2(Km+1:K)];
elseif Q2.K_sideNum<K
    beta_rad_2 = [beta_rad_2, beta_rad_1(Km+1:K)];
end


% 3layers sub trees
q_children_3Ls1 = cell(1,K);
q_children_3Ls2 = cell(1,K);

for k=1:Km
    q_children_3Ls1{k} = Q1.q_children{k};
    q_children_3Ls2{k} = Q2.q_children{k};
end

if Q1.K_sideNum<K
    for k=Km+1:K
        q_children_3Ls2{k} = Q2.q_children{k};
        q_children_3Ls1{k} = makeOneVirtualTree_rad_3layers(q_children_3Ls2{k});
    end
elseif Q2.K_sideNum<K
    for k=Km+1:K
        q_children_3Ls1{k} = Q1.q_children{k};
        q_children_3Ls2{k} = makeOneVirtualTree_rad_3layers(q_children_3Ls1{k});
    end
end

% make return values
Q1p = make_qCompTree_rad_4layers(q0_1, q_children_3Ls1, sk_1, Q1.b00_startP, ...
                                beta0_rad_1, beta_rad_1);
                            
Q2p = make_qCompTree_rad_4layers(q0_2, q_children_3Ls2, sk_2, Q2.b00_startP, ...
                                beta0_rad_2, beta_rad_2);

end