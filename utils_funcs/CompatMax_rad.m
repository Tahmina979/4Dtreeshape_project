function [Q1p,Q2p] = CompatMax_rad(Q1,Q2)


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


T = max( Q1.T_sidePointNums(1:Km), Q2.T_sidePointNums(1:Km) );
if Q1.K_sideNum>Km
    T = [T, Q1.T_sidePointNums(Km+1:end)];
elseif Q2.K_sideNum>Km
    T = [T, Q2.T_sidePointNums(Km+1:end)];
end

% beta rad
beta_rad_1 = Q1.beta_rad;
beta_rad_2 = Q2.beta_rad;
if Q1.K_sideNum<K
    beta_rad_1 = [beta_rad_1, beta_rad_2(Km+1:K)];
elseif Q2.K_sideNum<K
    beta_rad_2 = [beta_rad_2, beta_rad_1(Km+1:K)];
end


% q
q_1 = cell(1,K);
q_2 = cell(1,K);
for k=1:Km
%     T1 = Q1.T_sidePointNums(k); T2 = Q2.T_sidePointNums(k);
    T1 = size(Q1.q{k}, 2); T2 = size(Q2.q{k}, 2); maxT=max(T1, T2);
    tau = linspace(0,1,maxT);
    if T1 < maxT
        q_1{k} = interp1( linspace(0,1,T1),Q1.q{k}', tau )';
        beta_rad_1{k} = interp1(linspace(0,1,T1), Q1.beta_rad{k}, tau);
        q_2{k} = Q2.q{k};
    elseif T2 < maxT
        q_1{k} = Q1.q{k};
        q_2{k} = interp1( linspace(0,1,T2), Q2.q{k}', tau )';
        beta_rad_2{k} = interp1(linspace(0,1,T2), Q2.beta_rad{k}, tau);
    else
        q_1{k} = Q1.q{k};
        q_2{k} = Q2.q{k};
    end
end

if Q1.K_sideNum<K
    for k=Km+1:K
        Q1.dimension = 3;
        q_1{k} = zeros(Q1.dimension,size(Q2.q{k}, 2));
        q_2{k} = Q2.q{k};
    end
elseif Q2.K_sideNum<K
    for k=Km+1:K
        q_1{k} = Q1.q{k};
        Q2.dimension=3;
        q_2{k} = zeros(Q2.dimension,size(Q1.q{k}, 2));
    end
end

% sk
sk_1 = Q1.sk;
sk_2 = Q2.sk;
if Q1.K_sideNum<K
    sk_1 = [sk_1, sk_2(Km+1:K)];
    
elseif Q2.K_sideNum<K
    sk_2 = [sk_2, sk_1(Km+1:K)];
end

% make return values
Q1p = make_qST_rad(q0_1,q_1,sk_1,Q1.b00_startP, beta0_rad_1, beta_rad_1);
Q2p = make_qST_rad(q0_2,q_2,sk_2,Q2.b00_startP, beta0_rad_2, beta_rad_2);

end