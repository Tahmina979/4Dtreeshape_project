function [Mu,eigenVector, EVals,tree_pca,used_qCompTrees_Ready]=transfer_to_PCA_space_raw(ctree_1,ctree_2,ctree_3,ctree_4,ctree_5,ctree_6,ctree_7,num_within_seq1,num_within_seq2,num_within_seq3,num_within_seq4,num_within_seq5,num_within_seq6,num_within_seq7)
lam_m = 1; 
lam_s = 1;
lam_p = 1;

tNum=num_within_seq1+num_within_seq2+num_within_seq3+num_within_seq4+num_within_seq5+num_within_seq6+num_within_seq7;
used_qCompTrees=cell(1,tNum);

for i=1:num_within_seq1

  used_qCompTrees{i}=ctree_1{1,i};
end

p=num_within_seq1+1; p_end=num_within_seq1+num_within_seq2;

for i=p:p_end

  used_qCompTrees{i}=ctree_2{1,i-p+1};

end
p=p_end+1; p_end=p_end+num_within_seq3;
for i=p:p_end
    used_qCompTrees{i}=ctree_3{i-p+1};
end
p=p_end+1; p_end=p_end+num_within_seq4;
for i=p:p_end

    used_qCompTrees{i}=ctree_4{i-p+1};
end
p=p_end+1; p_end=p_end+num_within_seq5;
for i=p:p_end

    used_qCompTrees{i}=ctree_5{i-p+1};
end

p=p_end+1; p_end=p_end+num_within_seq6;
for i=p:p_end
  used_qCompTrees{i}=ctree_6{i-p+1};
end
p=p_end+1; p_end=p_end+num_within_seq7;
for i=p:p_end

    used_qCompTrees{i}=ctree_7{i-p+1};
end

[used_qCompTrees_Ready] = CompatMultiMax_rad_4layers(used_qCompTrees);


 
for i=1:tNum
   used_CompTrees_Ready{i}= qCompTree_to_CompTree_rad_4layers(used_qCompTrees_Ready{i});
 
end


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