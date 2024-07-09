function [ trunk1 ] = calcu_t_value_1_branch( trunk1 )

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

end

