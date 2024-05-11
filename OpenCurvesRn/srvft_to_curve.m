function p = srvft_to_curve(q)
%
% q - SRVFT rep of teh curve.
%   - q(:, 1:end-1) is the SRVF of the curve
%   - q(:, end) is its initial condition
%

qn      = q(:, 1:end-1);
[n,T]   = size(qn);
transVect = q(:, end);

for i = 1:T
    qnorm(i) = norm(qn(:,i));
end

for i = 1:n
    p(i,:) = [cumtrapz( qn(i,:).*qnorm )/(T) ] ;
end

p = p + repmat(transVect, 1, T);