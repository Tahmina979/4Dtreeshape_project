function [ST_1] = trunkRep_to_compTree_3layers_rad(trunk1)

%   
ST_1 = trunkRep_to_ST_rad(trunk1);

ST_1.beta_children = cell(1, numel(trunk1.children));

for i = 1: numel(trunk1.children)
    
    ST_1.beta_children{i} = trunkRep_to_ST_rad(trunk1.children{i});
    ST_1.beta_children{i}.beta_children = cell(1, numel(trunk1.children{i}.children));
    
    for j = 1: numel(ST_1.beta_children{i}.beta_children)
        ST_1.beta_children{i}.beta_children{j}.beta0 = ST_1.beta_children{i}.beta{j};
        
        Beta0 = ST_1.beta_children{i}.beta_children{j}.beta0;
        
        % --- Compute the function between the paras and radius ---
        t_paras = 0;
        total_sum = 0;     % --- calculate the total length of the truck.

        for jj=1: size(Beta0, 2)-1
            vec = Beta0(1:3, jj) - Beta0(1:3, jj+1);
            total_sum = total_sum + norm(vec, 2);
        end

        for ii=2: size(Beta0, 2)
            sum_length = 0;
            
            for jj=1: ii-1
                vec = Beta0(1:3, jj) - Beta0(1:3, jj+1);
                sum_length = sum_length + norm(vec, 2);
            end

            t_paras = [t_paras, sum_length/total_sum];
        end

        ST_1.beta_children{i}.beta_children{j}.t_paras = t_paras;
           
    end
    
    
    
end

end

