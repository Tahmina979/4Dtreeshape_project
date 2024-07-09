function [Q, locator_s, locator_e]= unflattenSimpleTree_rad(qX, lam_m, lam_s, lam_p, refQ, locator_s, locator_e)

% qX = zeros(n,p);
d = 3;

T0 = refQ.T0_pointNum;
dj = d*T0 +T0;
locator_e = locator_s + dj;

% store trunk and rad
q0_and_Rad = reshape(qX(locator_s:locator_e -1), [d+1, T0])*(T0-1);

q0 = q0_and_Rad(1:3, :)/sqrt(lam_m);
beta0_rad = q0_and_Rad(4, :);

% store subtree position
locator_s = locator_e;
locator_e = locator_s + numel(refQ.q);

sk = qX(locator_s:locator_e -1)/sqrt(lam_p);
sk(sk<0) = 0;

% make simple tree
q = cell(1, length(refQ.q));
Q = make_qST_rad(q0, q, sk,[0; 0 ;0], beta0_rad, refQ.beta_rad);

% unflatten subtrees
locator_s = locator_e;

for k=1: numel(refQ.q)

    T0 = size(refQ.q{k}, 2);
    dj = d*T0 +T0;
    locator_e = locator_s + dj;
    
    qk_and_rad = reshape(qX(locator_s:locator_e -1), [d+1, T0])*(T0-1);
    Q.q{k} = qk_and_rad(1:3, :)/sqrt(lam_s);
    Q.beta{k} = qk_and_rad(4, :);
    
    locator_s = locator_e;
    
end

            
end

