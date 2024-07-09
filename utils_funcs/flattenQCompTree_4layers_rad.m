function qX = flattenQCompTree_4layers_rad(qCompTree_4Ls, lam_m, lam_s, lam_p)

% qX = zeros(n,p);
d = 3;
locator_s= 0;
locator_e = 0;

locator_s = 1;
T0 = qCompTree_4Ls.T0_pointNum;
dj = d*T0 + T0;
locator_e = locator_s + dj;

% store trunk
Q0_and_Rad= [sqrt(lam_m)*qCompTree_4Ls.q0; qCompTree_4Ls.beta0_rad];
qX(locator_s:locator_e -1) = reshape(Q0_and_Rad, [1, dj])/(T0-1);

fprintf('%d-', locator_e-locator_s);

% store subtree position
locator_s = locator_e;
locator_e = locator_s + numel(qCompTree_4Ls.q_children);

qX(locator_s:locator_e -1) = sqrt(lam_p)*qCompTree_4Ls.sk(1:end);
fprintf('%d-', locator_e-locator_s);

% flatten subtrees
locator_s = locator_e;

for k=1: numel(qCompTree_4Ls.q_children)

    
    qX = flattenQCompTree_3layers_rad(qCompTree_4Ls.q_children{k}, lam_m, lam_s, lam_p, qX, locator_s);

    locator_s =length(qX)+1;
%     locator_e = locator_s+ numel(qX_children);
% 
%     qX(locator_s:locator_e-1) = qX_children;
%     
%     fprintf('%d-', locator_e-locator_s);
end

fprintf('\n');

            
end

