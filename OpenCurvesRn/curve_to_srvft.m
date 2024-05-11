% Converts the curve to its SRVFT representation, which is the SRVF and the
% initial condition so that the curve and its location (start point) can be reconstructed 
% 
% q   - the SRVFT: q(:, end-1) is the SRVF, q(:,end) is the initial condition
% len - the length of the curve
%
function [q,len] = curve_to_srvft(p, toNormalizeLength)

if nargin < 2
    toNormalizeLength = 1;
end

[n,N] = size(p);
for i = 1:n
    v(i,:) = gradient(p(i,:),1/(N));
end

% why dividing by N for the length????
len = sum(sqrt(sum(v.*v)))/N;

%% Normalize for scale (if needed)
if toNormalizeLength > 0
    v = v/len;        % This makes the length of the curve equal 1
end

%% Compute the SRVF
q = zeros(size(p, 1), size(p, 2) + 1);

for i = 1:N
    L = max(sqrt(norm(v(:,i))), 0.0001);
    q(:,i) = v(:,i)/L;
end

%% Append the initial condition, which results in the SRVFT rep
%q(:, end) = p(:, 1);