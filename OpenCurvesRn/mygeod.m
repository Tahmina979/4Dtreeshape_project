function [dist,X2n,q2n,X1,X2,q1, q2,gamI, distbefore]=mygeod(X1, X2, doRotation, doNormalizeScale)
%
% Distance is defined as the angle between the two curves. This is only
% valid for curves of length 1. Thus:
% - It is computed as acos when doNormalizeScale == 1, or Euclidean
% otheriwse.
%


%input: two curves X1 and X2 as Nxn (ND curve) 
% First dimension is the space dimension and second dimension is time
%   doRotation a flag indicatingwhether we also solve for optimal rotation
%
%output: 
% dist   - distance, also will display when you run the program; 
% X2n    - optimally registered curve X2; 
% q2n    - same as X2n except in q-fun form; 
% X1     - normalized curve X1; 
% q1     - same as X1 except in q-fun form;
% gamI   - Optimal diffo to apply to X2 in order to align it onto X1
%
%

if nargin < 3
    doRotation = 1;
end

if nargin < 4
    doNormalizeScale =0;
end

% Load some parameters, no need to change this
lam = 0;
    
% Resample the curves to have N points
N = size(X1, 2); %200;
    
% X1 = ReSampleCurve(X1,N);
if size(X2, 2) ~= N    
 % X2 = ReSampleCurve(X2,N); % done in 2nd phase of spatial registration
  
end

%Center curves, not really needed but good for display purposes
%     X1 = X1 - repmat(mean(X1')',1,size(X1,2));
%     X2 = X2 - repmat(mean(X2')',1,size(X2,2));    
    
%% Form the q function for representing curves and find best rotation
[q1] = curve_to_q(X1);%, doNormalizeScale);
[q2] = curve_to_q(X2);%, doNormalizeScale);
[n]  = size(q1,1);
    
%% Normalize for rotation (if needed)
if doRotation > 0   
   A = q1*q2';
   [U,S,V] = svd(A);
    if det(A) > 0
        S = eye(n);
    else
        S = eye(n);
        S(:,end) = -S(:,end);
    end
   Ot = U*S*V';
    %X2 = Ot*X2;
    %q2 = Ot*q2;
 end
 
% Applying optimal re-parameterization to the second curve
[gam] = DynamicProgrammingQ(q1/sqrt(InnerProd_Q(q1,q1)),q2/sqrt(InnerProd_Q(q2,q2)),lam,0);
gamI  = invertGamma(gam);
gamI  = (gamI-gamI(1))/(gamI(end)-gamI(1));
X2n   = Group_Action_by_Gamma_Coord(X2,gamI);
q2n   = curve_to_q(X2n);%,doNormalizeScale);

    
if doRotation > 0  
% Find optimal rotation
    A = q1*q2n';
    A = q1*q2';
    [U,S,V] = svd(A);
    if det(A) > 0
        S = eye(n);
    else
        S = eye(n);
        S(:,end) = -S(:,end);
    end
    Ot = U*S*V';
    %X2n = Ot*X2n;
    %q2n = Ot*q2n;
end

% Computing geodesic distance between the registered curves
N    = size(X1,2);
if doNormalizeScale == 1
    dist = acos(sum(sum(q1.*q2n))/N);   % arc length distance
    if nargout == 8
        distbefore = acos(sum(sum(q1.*q2))/N);
    end
else
        dist = sqrt(InnerProd_Q((q1-q2n),(q1-q2n)));  % sqrt(InnerProd_Q_closed(q1-q2n, q1-q2n)); 
   
   if nargout ==9
         distbefore = sqrt(InnerProd_Q(abs(q1-q2), abs(q1-q2))); % sqrt(InnerProd_Q_closed(q1-q2, q1-q2)); 
    end
end


 %sprintf('The distance between the two curves is %0.3f',dist);





