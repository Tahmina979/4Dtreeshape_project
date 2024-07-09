function [amount_layer,amount_branch] = Compute_treeInfo(compTree)

amount_layer=0;
amount_branch = 0;

amount_branch=amount_branch+1;  % the trunk

for i=1:numel(compTree.beta_children)
    amount_branch = amount_branch +1;
    for j=1:numel(compTree.beta_children{i}.beta_children)
        amount_branch = amount_branch +1;
        for k=1: numel(compTree.beta_children{i}.beta_children{j}.beta_children)
            amount_branch = amount_branch +1;
        end
    end
end

% layer amount

amount_layer4_branch =0;
for i=1:numel(compTree.beta_children)
    for j=1:numel(compTree.beta_children{i}.beta_children)
        amount_layer4_branch = amount_layer4_branch + ...
                            numel(compTree.beta_children{i}.beta_children{j}.beta);
        
    end
end

if amount_layer4_branch == 0
    amount_layer = 3;
else
    amount_layer = 4;
end






end