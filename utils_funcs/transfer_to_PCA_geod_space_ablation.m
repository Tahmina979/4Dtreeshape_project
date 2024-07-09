function [Mu,eigenVector, EVals,tree_pca,used_CompTrees_Ready,qX]=transfer_to_PCA_geod_space_ablation(aligned_tree_1,aligned_tree_2,num_within_seq1,num_within_seq2)
lam_m = 1; 
lam_s = 1;
lam_p = 1;

tNum=num_within_seq1+num_within_seq2;%+num_within_seq3+num_within_seq4+num_within_seq5+num_within_seq6+num_within_seq7;
used_qCompTrees=cell(1,tNum);

for i=1:num_within_seq1

  %used_qCompTrees{i}=CompTree_to_qCompTree_rad_4layers(aligned_tree_1{i});
  used_qCompTrees{i}=(aligned_tree_1{i});
end

p=num_within_seq1+1; p_end=tNum;

for i=p:p_end

  %used_qCompTrees{i}=CompTree_to_qCompTree_rad_4layers(aligned_tree_2{i-p+1});
  used_qCompTrees{i}=(aligned_tree_2{i-p+1});

end

[used_qCompTrees_Ready] = CompatMultiMax_rad_4layers(used_qCompTrees);


 

for i=1:tNum
    used_CompTrees_Ready{i}=qCompTree_to_CompTree_rad_4layers(used_qCompTrees_Ready{i});
end


qX = [];

for i = 1:tNum
 
   qX(i, :) = flattenphenoTree_4layers_rad(used_CompTrees_Ready{i},lam_p,lam_m,lam_s);
end


qXT= qX';

[Mu,eigenVector, EVals] = performEigenAnalysis(qXT);

%eigenVector=E(:,1:35);

for i=1:tNum

sample{i}=qXT(:,i)-Mu;

tree_pca(i,:)=sample{i}'* eigenVector;

end


end