% Input:
%   q is (n x Tq) matrix representing SRVF of curve in R^n
%   g is (1 x Tg) matrix representing reparameterization
% Output:
%   qg is (n x Tg) matrix, result of action (q,g) on SRVF
function qg = GammaActionQ(q,g)

[n,Tq] = size(q);
t = linspace(0,1,Tq);

dt = 1/(numel(g)-1);
gdot = gradient(g,dt);

qg = interp1(t,q',g)' .* repmat(sqrt(gdot),n,1);

end