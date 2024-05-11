function gamI = invertGamma(gam)

N = length(gam);
x = [1:N]/(N);

gamI = interp1(gam, x, x, 'linear');
