function compTree_mess = disorderSubtrees(compTree, rand_seed)

rng(rand_seed);
randIdxes = randperm(numel(compTree.beta));

% --- re-organzie compTree_mess ---
compTree_mess = compTree;


compTree_mess.T_sidePointNums = compTree.T_sidePointNums(randIdxes);
compTree_mess.tk_sideLocs = compTree.tk_sideLocs(randIdxes);
compTree_mess.beta = compTree.beta(randIdxes);
compTree_mess.beta_rad = compTree.beta_rad(randIdxes);
compTree_mess.beta_children = compTree.beta_children(randIdxes);

end





