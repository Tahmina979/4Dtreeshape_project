function Q = unflattenCompTree_4layers_rad(qX, lam_m, lam_s, lam_p, refQ)

% qX = zeros(n,p);
d = 3;

locator_s = 1;
T0 = refQ.T0_pointNum;
dj = d*T0 + T0;
locator_e = locator_s + dj;

% store trunk
q0_and_Rad = reshape(qX(locator_s:locator_e -1), [d+1, T0])*(T0-1);
q0 = q0_and_Rad(1:3, :)/sqrt(lam_m);
beta0_rad = q0_and_Rad(4, :);

% store subtree position
locator_s = locator_e;
locator_e = locator_s + numel(refQ.q_children);

sk = qX(locator_s:locator_e -1)/sqrt(lam_p);
sk(sk<0) = 0;

% Make 4 layers tree
% q_children_3Ls = cell(1, length(refQ.q_children));
q_children_3Ls = repmat({struct('q0', [])}, 1, length(refQ.q_children));


Q = make_qCompTree_rad_4layers(q0, q_children_3Ls, sk, [0; 0 ;0], ...
                                beta0_rad, refQ.beta_rad);
% unflatten subtrees
locator_s = locator_e;

for k=1: numel(Q.q_children)

    
    [Q.q_children{k}, locator_s,locator_e] = unflattenCompTree_3layers_rad(qX, lam_m, lam_s, lam_p, ...
                                                                           refQ.q_children{k}, locator_s, locator_e);
 
    Q.q{k} = Q.q_children{k}.q0;
    Q.beta{k} = Q.q_children{k}.beta0_rad;
end

            
end

