function [G,Q1p] = ReparamNoPerm_qST(Q1,Q2, lam_m,lam_s,lam_p)

% reparam main
[G.gam0,G.E] = DPQ_difflen(Q1.q0,Q2.q0);
q0 = GammaActionQ(Q1.q0,G.gam0);
G.E = lam_m*G.E;

% reparam first minK sides
minK = min(Q1.K,Q2.K);
K = max(Q1.K,Q2.K);
G.gam = cell(1,K);
q = cell(1,K);
for k=1:minK
    [G.gam{k},Eside] = DPQ_difflen(Q1.q{k},Q2.q{k});
    q{k} = GammaActionQ(Q1.q{k},G.gam{k});
    G.E = G.E + lam_s*Eside + lam_p*(Q1.sk(k)-Q2.sk(k))^2;
end
% arclen positions of first Q1.K
sk = Q1.sk;
% remaining qk and sk
if minK < Q2.K
    kkill = minK+1:Q2.K;
    for k=kkill
        G.gam{k} = linspace(0,1,Q2.T(k));
        q{k} = zeros(Q1.d,Q2.T(k));
    end
    sk = [sk, Q2.sk(kkill)];
    G.E = G.E + lam_s*sum(Q2.len(kkill));
elseif minK < Q1.K
    kkill = minK+1:Q1.K;
    for k=kkill
        G.gam{k} = linspace(0,1,Q1.T(k));
    end
    q(kkill) = Q1.q(kkill);
    G.E = G.E + lam_s*sum(Q1.len(kkill));
end

G.matched = [];
G.unmatched1 = 1:Q1.K;
G.unmatched2 = 1:Q2.K;

% [numel(q),numel(sk)]
Q1p = make_qST(q0,q,sk,Q1.b00);

