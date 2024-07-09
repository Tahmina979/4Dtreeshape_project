function [gam, E] = DPQ_difflen(q1,q2,t1,t2)

[d,T1] = size(q1);
[~,T2] = size(q2);


if nargin<3
    t1 = linspace(0,1,T1);
    t2 = linspace(0,1,T2);
end

% % put initial parameterization on a
% % common set of domain sample points
TT = min( 2*max(T1,T2)-1, lcm(T1-1,T2-1)+1 );
tt = linspace(0,1,TT);

% qq1 = interp1(t1,q1',tt)';
% qq2 = interp1(t2,q2',tt)';
qq1 = spline(t1,q1,tt);
qq2 = spline(t2,q2,tt);

ggam = DynamicProgrammingQ(qq1,qq2,0,0);

rtggdot = sqrt( gradient(ggam)./gradient(tt) );
% qq1gg = repmat(rtggdot,d,1).*interp1(tt,qq1', ggam)';
qq1gg = repmat(rtggdot,d,1).*spline(tt,qq1, ggam);
E = trapz(tt, sum((qq2-qq1gg).^2,1) );
% output is same length as t2,q2
gam = interp1(tt,ggam,t2);

end