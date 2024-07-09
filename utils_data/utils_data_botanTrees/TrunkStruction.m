function [trunk1] = TrunkStruction(branch1, branch_num1)

% get the trunk structure represention of a tree.

[branch1, bifurcation_num1] = calcu_branch_bifurcation(branch1, branch_num1);

trunk1 = branch1(1, 1);
bifur_num = length(branch1(1, 1).bifurcation);

trunk1.children = { };
for i = 1: bifur_num
   
    trunk1.children{i} = branch1(2, branch1(1, 1).bifurcation{i}.child);
    trunk1.children{i} = trunk1.children{i};  % �������ӵĻ���ȡ��һ��
      
    Num = numel(trunk1.children{i});
    for ii=1:Num   %same point of a branch in a layer could have two or more child
        trunk1.children{i}(ii).children = { };
        
        for j = 1: length(trunk1.children{i}(ii).bifurcation)
            trunk1.children{i}(ii).children{j} = branch1(3, trunk1.children{i}(ii).bifurcation{j}.child);   
        end
    end
    
end

% ���������ӵ� ��Ҫ���ӷֲ��ĸ���  �˴���������Bug
trunk1 = AdjustDataStruction(trunk1);

trunk1 = calcu_t_value(trunk1);

end

function [ branch, bifurcation_num ] = calcu_branch_bifurcation( branch, branch_num )
%CALCUBRANCHBIFUICATION Summary of this function goes here
%   Detailed explanation goes here
for i = 1: length(branch(1,1).point)
    branch(1, 1).point(i).child = [];
end

for i = 1: branch_num(2)
    for j = 1: length(branch(2,i).point)
        branch(2,i).point(j).child = [];
    end
end

for i = 1: branch_num(3)
    for j=1: length(branch(3,i).point)
        branch(3,i).point(j).child = [];
    end
end



% ������һ���֦�ɵ�ĺ�����������һ��ĵڼ���֦�ɵĵڼ����㡣 
for i = 1: branch_num(2)
    
    f_branch_id = branch(2, i).father_branch_id + 1;
    f_point_id = branch(2, i).father_point_id + 1;
    
    branch(1, f_branch_id).point(f_point_id).child = ...
                    [branch(1, f_branch_id).point(f_point_id).child, i];
end

% �����ڶ��������֦�ɵ�ĺ���������
% �ӵ�����֦�ɵĽǶȼ��õ��ڶ���֦�ɵĺ�����������
for i = 1: branch_num(3)
    
    f_branch_id = branch(3, i).father_branch_id + 1;
    f_point_id = branch(3, i).father_point_id + 1;
    
    branch(2, f_branch_id).point(f_point_id).child = ...
                    [branch(2, f_branch_id).point(f_point_id).child, i];
end

for i=1: branch_num(4)
    f_branch_id = branch(4, i).father_branch_id +1;
    f_point_id = branch(4, i).father_point_id +1;
    
    branch(3, f_branch_id).point(f_point_id).child = ...
                    [branch(3, f_branch_id).point(f_point_id).child, i];
end


% -------------------------------------------------------------------------

% 
 k = 1;
for i= 1: length(branch(1,1).point)
  
    % check the bifurcation point of the level
    if isempty (branch(1,1).point(i).child) ==0
        branch(1,1).bifurcation{k}.x = branch(1,1).point(i).x;
        branch(1,1).bifurcation{k}.y = branch(1,1).point(i).y;
        branch(1,1).bifurcation{k}.z = branch(1,1).point(i).z;
        branch(1,1).bifurcation{k}.r = branch(1,1).point(i).r;
        
        branch(1,1).bifurcation{k}.child = branch(1,1).point(i).child;
        
        branch(1,1).bifurcation{k}.origin_index = i;
        k = k + 1;
    end
    
end


% find the bifurcation points of the second layer branches
for i =1: branch_num(2)
    k = 1;
    for j = 1: length(branch(2,i).point)
        
        % check the bifurcation point of the branch
        if isempty (branch(2, i).point(j).child) ==0
            branch(2, i).bifurcation{k}.x = branch(2, i).point(j).x;
            branch(2, i).bifurcation{k}.y = branch(2, i).point(j).y;
            branch(2, i).bifurcation{k}.z = branch(2, i).point(j).z;
            branch(2, i).bifurcation{k}.r = branch(2, i).point(j).r;
            
            branch(2, i).bifurcation{k}.child = branch(2, i).point(j).child;

            branch(2, i).bifurcation{k}.origin_index = j;
            k = k + 1;
        end
    end
end

% find the bifurcation points of the third layer branches

for i=1: branch_num(3)
    k = 1;
    for j=1: length(branch(3, i).point)
         
        if isempty(branch(3, i).point(j).child) ==0
            branch(3, i).bifurcation{k}.x = branch(3, i).point(j).x;
            branch(3, i).bifurcation{k}.y = branch(3, i).point(j).y;
            branch(3, i).bifurcation{k}.z = branch(3, i).point(j).z;
            branch(3, i).bifurcation{k}.r = branch(3, i).point(j).r;
            
            branch(3, i).bifurcation{k}.child = branch(3, i).point(j).child;
            
            branch(3, i).bifurcation{k}.origin_index = j;
            k = k+1;
        end
    end
end
            

% %------------------------------------------------------------------------
for i=1: 3
    for j =1: branch_num(i)
        bifurcation_num(i, j) = length(branch(i, j).bifurcation);
    end
end


end


% here is the AdjustDataStruction function
function [ trunk1 ] = AdjustDataStruction( trunk1 )


add_bifur_num = 0;
k = 1;
new_children = { };
new_bifurcation = { };

bifur_num = numel(trunk1.bifurcation);
count = 1;
for i=1: bifur_num
%     if numel(trunk1.children{i}) == 1
%         new_children{k} = trunk1.children{i};
%         new_bifurcation{k} = trunk1.bifurcation{i};
%         k = k+1;
%     elseif numel(trunk1.children{i}) ==2
%         new_children{k} = trunk1.children{i}(1);
%         new_bifurcation{k} = trunk1.bifurcation{i}; k = k+1;
%         new_children{k} = trunk1.children{i}(2);
%         new_bifurcation{k} = trunk1.bifurcation{i}; k = k+1;
%     end
    
    for j=1: numel(trunk1.children{i})
        new_children{count} = trunk1.children{i}(j);
        new_bifurcation{count} = trunk1.bifurcation{i};
        count = count+1;
    end
    
end

trunk1.children = new_children;
trunk1.bifurcation = new_bifurcation;

for i=1: numel(trunk1.children)
    k = 1;
    new_children = { };
    new_bifurcation = { };
    for j = 1: numel(trunk1.children{i}.children)
        if numel(trunk1.children{i}.children{j}) == 1
            new_children{k} = trunk1.children{i}.children{j};
            new_bifurcation{k} = trunk1.children{i}.bifurcation{j};
            k = k+1;
        elseif numel(trunk1.children{i}.children{j}) ==2
            new_children{k} = trunk1.children{i}.children{j}(1);
            new_bifurcation{k} = trunk1.children{i}.bifurcation{j};
            k = k+1;
            new_children{k} = trunk1.children{i}.children{j}(2);
            new_bifurcation{k} = trunk1.children{i}.bifurcation{j};
            k = k+1;
            
        elseif numel(trunk1.children{i}.children{j}) ==3
            new_children{k} = trunk1.children{i}.children{j}(1);
            new_bifurcation{k} = trunk1.children{i}.bifurcation{j};
            k = k+1;
            new_children{k} = trunk1.children{i}.children{j}(2);
            new_bifurcation{k} = trunk1.children{i}.bifurcation{j};
            k = k+1;
            new_children{k} = trunk1.children{i}.children{j}(3);
            new_bifurcation{k} = trunk1.children{i}.bifurcation{j};
            k = k+1;
        end
        
    end
    
    trunk1.children{i}.children = new_children;
    trunk1.children{i}.bifurcation = new_bifurcation;
    
end
end


% ����ÿ���ֲ��Ĳ���tֵ
function [ trunk1 ] = calcu_t_value( trunk1 )
% calculate the t_values of the first layer

branch_data = [trunk1.point.x; trunk1.point.y; trunk1.point.z; trunk1.point.r];
SplineFunc= cscvn(branch_data);
length = SplineFunc.breaks(end);
for i=1: numel(trunk1.bifurcation)
    origin_ind = trunk1.bifurcation{i}.origin_index;
    trunk1.bifurcation{i}.t_value = SplineFunc.breaks(origin_ind)/ length;
end

trunk1.SplineFunc = SplineFunc;
trunk1.length = length;

%--------------------------------------------------------------------------
% calculate the t_value of the sencond layer
for i=1: numel(trunk1.children)
    
    clear branch_data SplineFunc;
    
    branch_data = [trunk1.children{i}.point.x; trunk1.children{i}.point.y; trunk1.children{i}.point.z; trunk1.children{i}.point.r];
    SplineFunc= cscvn(branch_data);
    length = SplineFunc.breaks(end);
    for j=1: numel(trunk1.children{i}.bifurcation)
        origin_ind = trunk1.children{i}.bifurcation{j}.origin_index;
        trunk1.children{i}.bifurcation{j}.t_value = SplineFunc.breaks(origin_ind)/ length;
    end
    trunk1.children{i}.SplineFunc = SplineFunc;
    trunk1.children{i}.length = length;
end

%% %%˳��Ҫ�ѵ����������������������ֵ�������

for i=1: numel(trunk1.children)
    for j=1: numel(trunk1.children{i}.children)
        clear branch_data SplineFunc;

        branch_data = [trunk1.children{i}.children{j}(1).point.x; ...
                            trunk1.children{i}.children{j}(1).point.y; ...
                                    trunk1.children{i}.children{j}(1).point.z; ...
                                            trunk1.children{i}.children{j}(1).point.r];
        SplineFunc= cscvn(branch_data);
        length = SplineFunc.breaks(end);

        trunk1.children{i}.children{j}(1).SplineFunc = SplineFunc;
        trunk1.children{i}.children{j}(1).length = length;
    end
end

end
