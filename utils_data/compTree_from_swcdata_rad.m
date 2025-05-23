function B = compTree_from_swcdata_rad(raw,ctype)

N = size(raw,1);

ap_id = -ones(N,2);

n = 1;
for i=1:N
    if raw(i,2)==ctype
        ap_id(i,1) = n;
        ap_id(i,2) = ap_id( raw(i,7), 1 );
        n = n+1;
    end
end


is_ap = (raw(:,2)==ctype);

raw_apical = raw( is_ap, : ); 



raw_apical(:,[1,7]) = ap_id(is_ap,:);

br=0;

raw_apical = raw;

first=0;
N = size(raw_apical,1);

% extract branches
beta = {};
diam = {};

branch_id = -ones(N,1);
branch_id(1) = 0;
prnt = [];
k=0;

% figure; 
% axis equal; hold on;
while any( branch_id == -1 )
    % get len to origin of each node
    len = zeros(N,1);
    for i=1:N
        p_i = raw_apical(i,7);
        if p_i == -1
            len(i) = 0;
        else
            xyz = raw_apical(i,3:5);
            p_xyz = raw_apical(p_i,3:5);
            len(i) = len(p_i) + norm(xyz-p_xyz);
        end
    end
   
    [~,end_id] = max(len);

  
%     disp(max(len));

    % extract next curve as longest path from origin or existing curve

    d = zeros(1,0);
    b = zeros(4,0);

    i = end_id;
    p_i = raw_apical(i,7);
    count=0;
    if first==0
        while count<48
         count=count+1;
        % i = p_i;
        % p_i = raw_apical(i,7);
        end
    end 
    while p_i ~= -1
        branch_id(i) = k;
        b = [raw_apical(i,3:6)', b];
        d = [raw_apical(i,6)',d];
        raw_apical(i,7) = -1;
        i = p_i;
        p_i = raw_apical(i,7); 
   
    end
    
   if first==0
    first=1;
    b = [raw_apical(i,3:6)', b];
    d = [raw_apical(i,6)',d];
   else    
    b = [raw_apical(i,3:6)', b];
    d = [raw_apical(i,6)',d];
   end

   if(size(b,2)==44)
   % b = flip(b,2);
    %d = flip(d,2);
   end
   
    if(size(b,2)<=4)
        %continue;
    end
   

    
    if k>0
        prnt(k) = branch_id(i);
    end
    k = k+1;
%     if k==20
%         break
%     end
   
   
    beta = [beta, b]; 
%     if isempty(beta) == 0
%         plot3(b(1, :), b(2, :), b(3, :), 'k', 'LineWidth', 4); hold on;
%     end
    diam = [diam, d];
   

end

prnt = [-1, prnt]+1;

lvl1 = find(prnt==1);

beta0 = beta{1};



[d,T0] = size(beta0);

t = linspace(0,1,T0);

% --- here just consider the second layers branches.
beta_2nd = beta(lvl1);

% --- draw beta_2nd ---
% for i=1: numel(beta_2nd)
%     plot3(beta_2nd{i}(1, :), beta_2nd{i}(2, :), beta_2nd{i}(3, :), 'r', 'LineWidth', 4);
%     hold on;
% end


K = numel(beta);

T = zeros(1,K);
tk = zeros(1,K);
for k=1:K
    T(k) = size(beta{k},2);
    bdist = sum( (beta0(1:3,:) -repmat(beta{k}(1:3,1),1,T0) ).^2, 1 );
    [~,i] = min(bdist);
    tk(k) = t(i);
    %break;
end


B = struct('t_paras',t, 'beta0',beta0, 'T0_pointNum',T0, 'K_sideNum',K,...
    'beta',{beta_2nd}, 'T_sidePointNums',T, 'tk_sideLocs',tk, 'dimension',d);



beta_id_left = find(prnt>1);

% --- here we construct a children index matrix ---
children_inds = cell(1, numel(beta));
for i=1: numel(beta)
    k = 1;
    for t = 1: numel(prnt)
        if prnt(t) == i
            children_inds{i}(k) = t;
            k = k+1;
        end
    end
end


%---  organize the tree structure ---

B = CompStruFieldValues_rad(B, beta, children_inds, 1);
% B.orig_ind = 1;
    
for i=1: numel(children_inds{1})
    ind = children_inds{1}(i);
    B.beta_children{i} = struct('t_paras',0, 'beta0',zeros(1, 0), 'T0_pointNum',0, 'K_sideNum',0,...
    'beta',cell(1,1), 'T_sidePointNums',0, 'tk_sideLocs',0, 'dimension',3, 'beta_children', cell(1,1));

    B.beta_children{i}= CompStruFieldValues_rad(B.beta_children{i}, beta, children_inds, ind);
    
    for t=1: numel(children_inds{ind})
        ind_ind = children_inds{ind}(t);
        B.beta_children{i}.beta_children{t} = struct('t_paras',0, 'beta0',zeros(1, 0), 'T0_pointNum',0, 'K_sideNum',0,...
        'beta',cell(1,1), 'T_sidePointNums',0, 'tk_sideLocs',0, 'dimension',3, 'beta_children', cell(1,1));
    
        B.beta_children{i}.beta_children{t}= CompStruFieldValues_rad(B.beta_children{i}.beta_children{t}, beta, children_inds, ind_ind);

        for k = 1: numel(children_inds{ind_ind})
            ind_ind_ind = children_inds{ind_ind}(k);
            
            B.beta_children{i}.beta_children{t}.beta_children{k} = struct('t_paras',0, 'beta0',zeros(1, 0), 'T0_pointNum',0, 'K_sideNum',0,...
            'beta',cell(1,1), 'T_sidePointNums',0, 'tk_sideLocs',0, 'dimension',3, 'beta_children', cell(1,1));
        
            B.beta_children{i}.beta_children{t}.beta_children{k}= ......
                    CompStruFieldValues_rad(B.beta_children{i}.beta_children{t}.beta_children{k}, beta, children_inds, ind_ind_ind);
        end
    end
    
end


% --- draw this representation ---
% --- first generate the colors ---
c = distinguishable_colors(5,[.8,.8,.8;0,0,0],[20,20]);




% separate radius
B = SepRadiusData(B);

for i = 1: B.K_sideNum
    B.beta_children{i} = SepRadiusData(B.beta_children{i});
    
    for j=1:B.beta_children{i}.K_sideNum
        B.beta_children{i}.beta_children{j} = SepRadiusData(B.beta_children{i}.beta_children{j});
    end
    
end


end


function ST = SepRadiusData(ST)


ST.beta0_rad = zeros(1, 0);
ST.beta_rad = cell(1, ST.K_sideNum);

ST.beta0_rad = ST.beta0(4, :);
ST.beta0(4, :) = [];

for i = 1: ST.K_sideNum
    ST.beta_rad{i} = ST.beta{i}(4, :);
    ST.beta{i}(4, :) = [];
end


end

