function [Q1p,Q2p] = CompatMax(Q1,Q2)


%%% main branch
T0_pointNum = max(Q1.T0_pointNum, Q2.T0_pointNum);
t = linspace(0,1,T0_pointNum);

% q0
if Q1.T0_pointNum<T0_pointNum
    q0_1 = interp1(Q1.t_paras,Q1.q0', t)';
    q0_2 = Q2.q0;
elseif Q2.T0_pointNum<T0_pointNum
    q0_1 = Q1.q0;
    q0_2 = interp1(Q2.t_paras,Q2.q0', t)';
else
    q0_1 = Q1.q0;
    q0_2 = Q2.q0;
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

% q
q_1 = cell(1,K);
q_2 = cell(1,K);
for k=1:Km
%     T1 = Q1.T_sidePointNums(k); T2 = Q2.T_sidePointNums(k);
    T1 = length(Q1.q{k}); T2 = length(Q2.q{k});
    tau = linspace(0,1,T(k));
    if T1 < T(k)
        q_1{k} = interp1( linspace(0,1,T1),Q1.q{k}', tau )';
        q_2{k} = Q2.q{k};
    elseif T2 < T(k)
        q_1{k} = Q1.q{k};
        q_2{k} = interp1( linspace(0,1,T2),Q2.q{k}', tau )';
    else
        q_1{k} = Q1.q{k};
        q_2{k} = Q2.q{k};
    end
end

if Q1.K_sideNum<K
    for k=Km+1:K
        q_1{k} = zeros(Q1.dimension,T(k));
        q_2{k} = Q2.q{k};
    end
elseif Q2.K_sideNum<K
    for k=Km+1:K
        q_1{k} = Q1.q{k};
        q_2{k} = zeros(Q2.dimension,T(k));
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
Q1p = make_qST(q0_1,q_1,sk_1,Q1.b00_startP);
Q2p = make_qST(q0_2,q_2,sk_2,Q2.b00_startP);

end