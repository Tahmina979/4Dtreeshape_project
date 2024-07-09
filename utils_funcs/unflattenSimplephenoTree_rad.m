function [Q, locator_s, locator_e]= unflattenSimplephenoTree_rad(qX,lam_m, lam_s, lam_p, refQ, locator_s, locator_e)

% qX = zeros(n,p);
d = 3;

T0 = refQ.T0_pointNum;
dj = d*T0 +T0;
locator_e = locator_s + dj;

% store trunk and rad



% store trunk
beta0_and_Rad= reshape(qX(locator_s:locator_e -1), [d+1, T0])*(T0-1);
beta0 = beta0_and_Rad(1:3, :)/sqrt(lam_m);
beta0_rad = beta0_and_Rad(4, :);

% store subtree position
locator_s = locator_e;
locator_e = locator_s + numel(refQ.beta);

tk_sideLocs = qX(locator_s:locator_e -1)/sqrt(lam_p);
%tk_sideLocs(tk_sideLocs<0) = 0;

% make simple tree
beta = cell(1, length(refQ.beta));
%beta = repmat({struct('beta0', [])}, 1, length(refQ.beta));

Q = make_phenoST_rad(beta0, beta, tk_sideLocs,[0; 0 ;0], beta0_rad, refQ.beta_rad);

% unflatten subtrees
locator_s = locator_e;

for k=1: numel(refQ.beta)

    T0 = size(refQ.beta{k}, 2);
    dj = d*T0+T0;
    locator_e = locator_s + dj;
   
    betak_and_rad = reshape(qX(locator_s:locator_e -1), [d+1, T0])*(T0-1);
    Q.beta{k} = betak_and_rad(1:3, :)/sqrt(lam_s);
    Q.beta_rad{k} = betak_and_rad(4, :);
    
    locator_s = locator_e;
    
end
            
end

