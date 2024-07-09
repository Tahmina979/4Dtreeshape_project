function [Q1, Q2] = ResetqCompTree_q_and_beta_rad(Q1, Q2)

for i= 1: length(Q1.q_children)
    Q1.q{i} = Q1.q_children{i}.q0;
    Q1.beta_rad{i} = Q1.q_children{i}.beta0_rad; 
end

for i= 1: length(Q2.q_children)
    Q2.q{i} = Q2.q_children{i}.q0;
    Q2.beta_rad{i} = Q2.q_children{i}.beta0_rad; 
end


end