% Converts the curve to its SRVF representation
%
%
function [q] = curve_to_qo(p, toNormalizeLength)

if nargin < 2
    toNormalizeLength = 1;
end

[n,N] = size(p);
for i = 1:n
    v(i,:) = gradient(p(i,:),1/(N-1));
end

%len = sum(sqrt(sum(v.*v)))/N;

%% Normalize for scale (if needed)
if toNormalizeLength > 0
   % v   = v/len;        % This makes the length of the curve equal 1
end

%% Compute the SRVF
q = zeros(size(p));

for i = 1:N
    L(i) = max(sqrt(norm(v(:,i))), 0.0001);
    q(:,i) = v(:,i)/L(i);
end