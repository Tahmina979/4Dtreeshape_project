function [Mu, eigenVectors, EVals] = performEigenAnalysis(V, nModes, Mu)
%
% V is dim x N, where
% - dim: the dimensionality of the data
% - N  : number of samples
% - Mu: the mean - if not provided, it will be computed inside this
% function

[dim n] = size(V);

if nargin < 2,
    nModes = n;
end

% Mean shape
if nargin < 3,
    Mu = sum(V, 2) / n;
end
U = V - repmat(Mu, 1, n);

% Modes
if n<dim
    S = U' * U;
    [V D] = eigs(S,nModes-1);
   
    eigenVectors = U * V;
    
    % normalizing the eigenvectors
    for i=1:size(eigenVectors, 2),
        eigenVectors(:, i) = eigenVectors(:, i)  / norm(eigenVectors(:, i) );
    end
    EVals  = diag( D) / (n - 1); 
   
else
    % Covariance matrix
    C = U * U' / (n-1);
    % SVD decomposition
    [eigenVectors, D] = eigs(C, nModes); %, 50);
    EVals  = diag( D); 
    
end

% %% make the sign of the eigenvects consistent (i.e align them with axis dirs)
% for i=1:size(eigenVectors, 2),
%     w = zeros(size(eigenVectors, 1), 1);
%     w(i) = 1;   
% 
% %     disp(dot(eigenVectors(:, i), w) )
% %     pause
%     eigenVectors(:, i) = sign(dot(eigenVectors(:, i), w)) * eigenVectors(:, i)  / norm(eigenVectors(:, i) );
% end

