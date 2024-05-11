function val = InnerProd_Q(q1, q2)

[n, T] = size(q1);
val    = trapz(linspace(0,1,T), sum(q1.*q2));

return;