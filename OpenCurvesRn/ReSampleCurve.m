function Xn = ReSampleCurve(X,N)

[n,T] = size(X)
del(1) = 0;
for r = 2:T
    del(r) = norm(X(:,r) - X(:,r-1));
end
cumdel = cumsum(del)/sum(del);   

newdel = [0:N-1]/(N-1);

for j=1:n
    Xn(j,:) = interp1(cumdel,X(j,1:T),newdel,'linear');   
end

end