function [ST_1] = trunkRep_to_ST_rad(trunk1)

% 
for i = 1: numel(trunk1.point)
    ST_1.beta0(1, i) = trunk1.point(i).x;
    ST_1.beta0(2, i) = trunk1.point(i).y;
    ST_1.beta0(3, i) = trunk1.point(i).z;
    ST_1.beta0_rad(i) = trunk1.point(i).r;
end

ST_1.T0_pointNum = numel(trunk1.point);
ST_1.K_sideNum = numel(trunk1.children);

ST_1.T_sidePointNums = [];
ST_1.tk_sideLocs = [];
ST_1.beta = {};


ST_1.beta_rad={};
for i=1: numel(trunk1.children)
    for j=1: numel(trunk1.children{i}.point)
    ST_1.beta{i}(1, j) = trunk1.children{i}.point(j).x;
    ST_1.beta{i}(2, j) = trunk1.children{i}.point(j).y;
    ST_1.beta{i}(3, j) = trunk1.children{i}.point(j).z;
    ST_1.beta_rad{i}(j) = trunk1.children{i}.point(j).r;
    end
    
    ST_1.T_sidePointNums = [ST_1.T_sidePointNums, numel(trunk1.children{i}.point)];
    
    % --- compute the tk_sideLocs ---
    
%     ST_1.tk_sideLocs = [ST_1.tk_sideLocs, trunk1.bifurcation{i}.t_value];
end

ST_1.dimension = 3;
ST_1.t_paras = 0;
total_sum = 0;     % --- calculate the total length of the truck.

for j=1: size(ST_1.beta0, 2)-1
    vec = ST_1.beta0(1:3, j) - ST_1.beta0(1:3, j+1);
    total_sum = total_sum + norm(vec, 2);
end

for i=2: size(ST_1.beta0, 2)

    Mat1= ST_1.beta0(:, 1:i-1);
    Mat2= ST_1.beta0(:, 2:i);
    sum_length = sum(vecnorm(Mat1- Mat2) );
    
    ST_1.t_paras = [ST_1.t_paras, sum_length/total_sum];
end

% --- Compute the tk_sideLocs ---

for i=1: numel(trunk1.children)
    
    % --- compute the tk_sideLocs ---
     cur_bif = [trunk1.bifurcation{i}.x, trunk1.bifurcation{i}.y, trunk1.bifurcation{i}.z];
     
     k = 1;
     for j=1: size(ST_1.beta0, 2)
         
         if norm([ST_1.beta0(1: 3, j) - cur_bif']) < 0.000001
             k = j;
         end
     end
     
    ST_1.tk_sideLocs = [ST_1.tk_sideLocs, ST_1.t_paras(k)];
end


end

