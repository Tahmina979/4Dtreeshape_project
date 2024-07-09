% Algorithm 2 in paper
function [G,qST1p, qST2p] = ReparamPerm_qST_complex_forZero(qST1,qST2, lam_m,lam_s,lam_p)
% function [matched,Etotal,E] = ReparamPerm_qST(qST1,qST2, lam_m,lam_s,lam_p)

%%%%% MATCH SIDES %%%%% 
K1 = qST1.K;
K2 = qST2.K;

s1k = qST1.sk;
s2k = qST2.sk;

% build matching energy and side-branch reparameterization matrices


% % --- case 0 : parallel computing of the distance matrix ---
% AA0 = tic;
% % A2 = tic;
% gam_side_all = cell(K1,K2);
% E1 = zeros(K1,K2);
% for i=1:K1
%     %fprintf('%d...',i);
%     for j=1:K2
%         [gam_side_all{i,j}, Eside0] = DPQ_difflen(qST1.q{i},qST2.q{j});
%         [G, ~] = ReparamPerm_qST_new(qST1.q_children{i},qST2.q_children{j}, lam_m,lam_s,lam_p);
%         Eside = sqrt(G.E);
%         E0(i,j) = lam_s*Eside + lam_p*(s1k(i)-s2k(j)).^2;
%     end
%     
% %     q_zero = wg_root_makeZeroST(qST1.q_children{i});
% %     [G1] = ReparamPerm_qST_forZero(qST1.q_children{i}, q_zero, lam_m,lam_s,lam_p);
% %     E(i,K2+1:end) = lam_s* sqrt(G1.E);
% end
% 
% % for j=1:K2
% %     clear q_zero
% %     q_zero = wg_root_makeZeroST(qST2.q_children{j});
% %     [G2] = ReparamPerm_qST_forZero(q_zero, qST2.q_children{j}, lam_m,lam_s,lam_p);
% %     
% %     E(K1+1:end,j) = lam_s* sqrt(G2.E);
% % end
% 
% % t1 = toc(AA);
% 
% [Mvec0,Efin0] = munkres(E0);
% tt0 = toc(AA0);
% global tt0_t0;
% 
% % tt1_t1 = [tt1, t1];
% tt0_t0 = [tt0_t0,tt0];
% 

%% --- Original Method ---

gam_side_all = cell(K1,K2);

E = zeros(K1+K2,K1+K2);

for i=1:K1
    %fprintf('%d...',i);
    for j=1:K2
        % --- the process here just makes code safe. ---
        if isempty(qST1.q{i})
            qST1.q{i} = [0,0;0,0;0,0];
        end
        
        if isempty(qST2.q{j})
            qST2.q{j}= [0,0;0,0;0,0];
        end
        
        [gam_side_all{i,j}, Eside0] = DPQ_difflen(qST1.q{i}, qST2.q{j});
        [G, ~] = ReparamPerm_qST(qST1.q_children{i}, qST2.q_children{j}, lam_m,lam_s,lam_p);
        Eside = sqrt(G.E);
        
        E(i,j) = lam_s*Eside + lam_p*(s1k(i)-s2k(j)).^2;
    end
    
    q_zero = wg_root_makeZeroST(qST1.q_children{i});
    [G1] = ReparamPerm_qST_forZero(qST1.q_children{i}, q_zero, lam_m,lam_s,lam_p);
    E(i,K2+1:end) = lam_s* sqrt(G1.E);
end

for j=1:K2
    clear q_zero
    q_zero = wg_root_makeZeroST(qST2.q_children{j});
    [G2] = ReparamPerm_qST_forZero(q_zero, qST2.q_children{j}, lam_m,lam_s,lam_p);
    
    E(K1+1:end,j) = lam_s* sqrt(G2.E);
end
% t2 = toc(BB);

% optimize assignment
% a = tic;
[Mvec,Efin] = munkres(E);

% get list of matched pairs and unmatched from munkres output
matched = zeros(0,2);
unmatched1 = [];
unmatched2 = [];
gam_side = cell(1,0);

for i=1:K1
    if Mvec(i)<=K2
        matched = [matched; i,Mvec(i)];
        gam_side = [gam_side, gam_side_all{i,Mvec(i)}];
    else
        unmatched1 = [unmatched1, i];
    end
end

for i= K1+(1:K2)
    if Mvec(i)<=K2
        unmatched2 = [unmatched2, Mvec(i)];
    end
end


%% --- Improved Method ---

% gam_side_all = cell(K1,K2);
% E1 = zeros(K1,K2);
% for i=1:K1
%     %fprintf('%d...',i);
%     for j=1:K2
%         [gam_side_all{i,j}, Eside0] = DPQ_difflen(qST1.q{i},qST2.q{j});
%         [G, ~] = ReparamPerm_qST_new(qST1.q_children{i}, qST2.q_children{j}, lam_m,lam_s,lam_p);
%         Eside = sqrt(G.E);
%         E1(i,j) = lam_s*Eside + lam_p*(s1k(i)-s2k(j)).^2;
%     end
%     
% end
% 
% % t1 = toc(AA);
% 
% [Mvec1,Efin] = munkres(E1);
% 
% matched = zeros(0,2);
% unmatched1 = [];
% unmatched2 = [];
% gam_side = cell(1,0);
% 
% for i=1:K1
%     if Mvec1(i) ~=0
%         matched = [matched; i, Mvec1(i)];
%         gam_side = [gam_side, gam_side_all{i,Mvec1(i)}];
%     else
%         unmatched1 = [unmatched1, i];
%     end
% end
% 
% for i=1:K2
%     if any(Mvec1 == i) == 0
%         unmatched2 = [unmatched2, i];
%     end
% end

%%

%fprintf('\nAlign main...\n');

[gam0, Emain] = DPQ_difflen(qST1.q0, qST2.q0);  % --- DPQ_difflen computes the difference between two q-space branch.

Etotal = lam_m*Emain + Efin;

G = struct('E',Etotal, 'gam0',gam0, 'gam',{gam_side}, ...
    'matched',matched, 'unmatched1',unmatched1, 'unmatched2',unmatched2);

%%%%% TRANSFORM qST1 %%%%%
%%% q0, t



end

