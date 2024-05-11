function [dist,X2n,X1,X2,gamI, distbefore]=mygeod_org(X1, X2, doRotation, doNormalizeScale)
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
 % X2 = ReSampleCurve(X2,N); % not required, it is done for 2nd phase of spatial registration
end

%Center curves, not really needed but good for display purposes
%     X1 = X1 - repmat(mean(X1')',1,size(X1,2));
%     X2 = X2 - repmat(mean(X2')',1,size(X2,2));    

% Applying optimal re-parameterization to the second curve
[gam] = DynamicProgrammingQ(X1/sqrt(InnerProd_Q(X1,X1)),X2/sqrt(InnerProd_Q(X2,X2)),lam,0);
gamI  = invertGamma(gam);
gamI  = (gamI-gamI(1))/(gamI(end)-gamI(1));
X2n   = Group_Action_by_Gamma_Coord(X2,gamI);

    


% Computing geodesic distance between the registered curves
N    = size(X1,2);
if doNormalizeScale == 0
    dist = sqrt(InnerProd_Q((X1-X2n),(X1-X2n)));
end
        
   
   if nargout ==6
         distbefore = sqrt(InnerProd_Q(abs(X1-X2), abs(X1-X2))); % sqrt(InnerProd_Q_closed(q1-q2, q1-q2)); 
   end

end


 %sprintf('The distance between the two curves is %0.3f',dist);





