function [Q1p,Q2p] = CompatMax_rad_3layers_within(Q1, Q2)

ret_value=0;
%%% main branch
T0_pointNum = max(Q1.T0_pointNum, Q2.T0_pointNum);
T0_pointNum_min = min(Q1.T0_pointNum, Q2.T0_pointNum);
t = linspace(0,1,T0_pointNum);
t_min = linspace(0,1,T0_pointNum_min);

% q0
if Q1.T0_pointNum<T0_pointNum
    q0_1 = interp1(Q1.t_paras,Q1.q0', t)'; 
    beta0_rad_1 = interp1(Q1.t_paras, Q1.beta0_rad, t);
    
    q0_2 = Q2.q0;
    beta0_rad_2 = Q2.beta0_rad;
elseif Q2.T0_pointNum<T0_pointNum
    q0_1 = interp1(Q1.t_paras,Q1.q0', t_min)'; 
    beta0_rad_1 = interp1(Q1.t_paras, Q1.beta0_rad, t_min);
    
    q0_2 = Q2.q0;
    beta0_rad_2 = Q2.beta0_rad;
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
    ret_value=1;
end
if ret_value~=1
% beta rad
beta_rad_1 = Q1.beta_rad;
beta_rad_2 = Q2.beta_rad;

    beta_rad_1 = [beta_rad_1, beta_rad_2(Km+1:K)];


% 2layers sub trees
q_children_ST1 = cell(1,K);
q_children_ST2 = cell(1,K);

for k=1:Km
    q_children_ST1{k} = Q1.q_children{k};
    q_children_ST2{k} = Q2.q_children{k};
end


    for k=Km+1:K
        q_children_ST2{k} = Q2.q_children{k};
        q_children_ST1{k} = makeZeroST_rad(q_children_ST2{k});
    end


% make return values
Q1p = make_qCompTree_rad_3layers(q0_1, q_children_ST1, sk_1, Q1.b00_startP, ...
                                 beta0_rad_1, beta_rad_1);
                            
Q2p = make_qCompTree_rad_3layers(q0_2, q_children_ST2, sk_2, Q2.b00_startP, ...
                                 beta0_rad_2, beta_rad_2);
else
    Q1p=Q1;
    Q2p=Q2;
end

end