function [Mu,eigenVector, EVals,tree_pca,used_qCompTrees_Ready,qX]=transfer_to_PCA_temp_error_space(aligned_q_tree_1,aligned_q_tree_2,GT,num_within_seq1,num_within_seq2,GTnum)
lam_m = 1; 
lam_s = 1;
lam_p = 1;

tNum=num_within_seq1+num_within_seq2+GTnum;%+num_within_seq3+num_within_seq4+num_within_seq5+num_within_seq6+num_within_seq7;
used_qCompTrees=cell(1,tNum);

for i=1:num_within_seq1

  used_qCompTrees{i}=aligned_q_tree_1{i};
end

p=num_within_seq1+1; p_end=num_within_seq1+num_within_seq2;

for i=p:p_end

  used_qCompTrees{i}=aligned_q_tree_2{i-p+1};

end

p=p_end+1; p_end=tNum;

for i=p:p_end

  used_qCompTrees{i}=GT{i-p+1};

end

[used_qCompTrees_Ready] = CompatMultiMax_rad_4layers(used_qCompTrees);


 



qX = [];

for i = 1:tNum
 
   qX(i, :) = flattenQCompTree_4layers_rad(used_qCompTrees_Ready{i}, lam_m, lam_s, lam_p);
end


qXT= qX';

[Mu,eigenVector, EVals] = performEigenAnalysis(qXT);

%eigenVector=E(:,1:35);

for i=1:tNum

sample{i}=qXT(:,i)-Mu;

tree_pca(i,:)=sample{i}'* eigenVector;

end


end