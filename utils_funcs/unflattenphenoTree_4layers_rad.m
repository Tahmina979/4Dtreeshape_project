function Q = unflattenphenoTree_4layers_rad(qX,lam_m, lam_s, lam_p, refQ)

% qX = zeros(n,p);
d = 3;

locator_s = 1;
T0 = refQ.T0_pointNum;
dj = d*T0 + T0;
locator_e = locator_s + dj;

% store trunk

beta0_and_Rad = reshape(qX(locator_s:locator_e -1), [d+1, T0])*(T0-1);
beta0 = beta0_and_Rad(1:3, :)/sqrt(lam_m);
beta0_rad = beta0_and_Rad(4, :);

% store subtree position
locator_s = locator_e;
locator_e = locator_s + numel(refQ.beta_children);



tk_sideLocs = qX(locator_s:locator_e -1)/sqrt(lam_p);

%tk_sideLocs(tk_sideLocs<0) = 0;

% Make 4 layers tree
% q_children_3Ls = cell(1, length(refQ.q_children));
beta_children_3Ls = repmat({struct('beta0', [])}, 1, length(refQ.beta_children));


Q = make_phenoTree_rad_4layers(beta0, beta_children_3Ls, tk_sideLocs, [0; 0 ;0], ...
                                beta0_rad, refQ.beta_rad);
% unflatten subtrees
locator_s = locator_e;

for k=1: numel(Q.beta_children)

    
    [Q.beta_children{k}, locator_s,locator_e] = unflattenphenoTree_3layers_rad(qX, lam_m, lam_s, lam_p, ...
                                                                           refQ.beta_children{k}, locator_s, locator_e);
   
    Q.beta{k} = Q.beta_children{k}.beta0;
    Q.beta_rad{k} = Q.beta_children{k}.beta0_rad;
end

            
end

