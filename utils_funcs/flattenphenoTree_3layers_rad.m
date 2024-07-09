function qX = flattenphenoTree_3layers_rad(CompTree_3Ls, lam_m, lam_s, lam_p, qX, locator_s)


% qX = zeros(n,p);
d = 3;

locator_e = 0;

T0 = CompTree_3Ls.T0_pointNum;
dj = d*T0 + T0;
locator_e = locator_s + dj;

% store trunk
beta0_and_Rad= [sqrt(lam_m)*CompTree_3Ls.beta0; CompTree_3Ls.beta0_rad];
qX(locator_s:locator_e -1) = reshape(beta0_and_Rad, [1, dj])/(T0-1);
fprintf('3Ls%d-', locator_e-locator_s);

% store subtree position
locator_s = locator_e;
locator_e = locator_s + numel(CompTree_3Ls.beta_children);

qX(locator_s:locator_e -1) = sqrt(lam_p)*CompTree_3Ls.sk(1:end);
fprintf('3Ls%d-', locator_e-locator_s);

% flatten subtrees
locator_s = locator_e;

for k=1: numel(CompTree_3Ls.beta_children)

    qX = flattenphenoSimpleTree_rad(CompTree_3Ls.beta_children{k}, lam_m, lam_s, lam_p, qX, locator_s);
    locator_s =length(qX)+1;

%     locator_s =locator_e;
%     locator_e = locator_s+ numel(qX_children);

%     qX(locator_s:locator_e-1) = qX_children;
%     fprintf('Simp%d-%d-', k,locator_e-locator_s);
end








end