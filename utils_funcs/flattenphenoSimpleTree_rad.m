function qX = flattenphenoSimpleTree_rad(SimpleTree, lam_m, lam_s, lam_p, qX, locator_s)


% qX = zeros(n,p);
d = 3;
locator_e = 0;


T0 = SimpleTree.T0_pointNum;
dj = d*T0+ T0;
locator_e = locator_s + dj;

% store trunk
beta0_and_Rad= [sqrt(lam_m)*SimpleTree.beta0; SimpleTree.beta0_rad];
qX(locator_s:locator_e -1) = reshape(beta0_and_Rad, [1, dj])/(T0-1);

% store subtree position
locator_s = locator_e;
locator_e = locator_s + numel(SimpleTree.beta);

qX(locator_s:locator_e -1) = sqrt(lam_p)*SimpleTree.tk_sideLocs(1:end);


% flatten subtrees
locator_s = locator_e;

for k=1: numel(SimpleTree.beta)

    T0 = size(SimpleTree.beta{k}, 2);
    dj = d*T0 +T0;
    locator_e = locator_s + dj;

    % store trunk
    betak_and_Rad= [sqrt(lam_s)*SimpleTree.beta{k}; SimpleTree.beta_rad{k}];
    qX(locator_s:locator_e -1) = reshape(betak_and_Rad, [1, dj])/(T0-1);
    
    locator_s = length(qX)+1;
    

end




end