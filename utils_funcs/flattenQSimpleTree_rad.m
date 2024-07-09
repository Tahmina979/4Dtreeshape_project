function qX = flattenSimpleTree_rad(qSimpleTree, lam_m, lam_s, lam_p, qX, locator_s)


% qX = zeros(n,p);
d = 3;
locator_e = 0;


T0 = qSimpleTree.T0_pointNum;
dj = d*T0 + T0;
locator_e = locator_s + dj;

% store trunk
Q0_and_Rad= [sqrt(lam_m)*qSimpleTree.q0; qSimpleTree.beta0_rad];
qX(locator_s:locator_e -1) = reshape(Q0_and_Rad, [1, dj])/(T0-1);

% store subtree position
locator_s = locator_e;
locator_e = locator_s + numel(qSimpleTree.q);

qX(locator_s:locator_e -1) = sqrt(lam_p)*qSimpleTree.sk(1:end);


% flatten subtrees
locator_s = locator_e;

for k=1: numel(qSimpleTree.q)

    T0 = size(qSimpleTree.q{k}, 2);
    dj = d*T0 +T0;
    locator_e = locator_s + dj;

    % store trunk
    Qk_and_Rad= [sqrt(lam_s)*qSimpleTree.q{k}; qSimpleTree.beta_rad{k}];
    qX(locator_s:locator_e -1) = reshape(Qk_and_Rad, [1, dj])/(T0-1);
    
    locator_s = length(qX)+1;
    

end




end