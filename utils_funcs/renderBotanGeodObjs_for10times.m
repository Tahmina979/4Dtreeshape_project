clc; 
% clear; 
close all;
addpath('func_render/');
addpath('func_other/');
% mesh_dir = 'data/';
mesh_dir = [obj_folder, '/'];

% load the shape
N=9;
S = cell(1, N);
for i = 1:N
    S{i} = MESH_IO.read_shape([mesh_dir, num2str(i),'.obj']);
end

%% Draw geodesic
 renderOptions = {'RotationOps',{[-90,0,0],[0,0,0]},...  % Botanical trees: [-90,0,0],[0,0,0]
                  'CameraPos',[-0.1,10],...             % default: [-10, 10]
                  'FaceAlpha',0.9,...
                  'BackgroundColor',[0.9, 0.9, 0.9]}; % you can change the background color here
              
fig = figure('visible','on');
axis equal; hold on;
default_c = get(gca,'colororder');

for i=1:N

    M = S{i};
    color_botanTree = [0.2, 0.1, 0];
%     color_neuroTree = [0.0, 0.0, 0.5];
%     color_neuroTree = default_c(1, :);
    
    M_col = repmat(color_botanTree, M.nv, 1);    

    moveVec = [(i-1)*5.5, 0, 0];  % --- BotanTrees: [(i-1)*4, 0, 0], NeuroTrees; [(i-1)*4, 0, 0]
    
    tm1 = tic();
    [~,~, S1_new] = render_mesh(M,'MeshVtxColor',M_col,...
                               'VtxPos',M.surface.VERT + repmat(moveVec, M.nv,1),... % translate the second shape such that S1 and S2 are not overlapped
                                renderOptions{:});
    hold on;
%     all_edges = get_edge_list(S{i});  % get the edge list
%     vertex = S1_new.surface.VERT';    % get the vertex positions after the rotations!
% 
%     
%     t = trimesh(S1_new.surface.TRIV, vertex(1,:), vertex(2,:), vertex(3,:),...
%                 'FaceColor','none',...
%                 'EdgeColor', color_botanTree,...
%                 'LineWidth', 0.2);
% 
%     hold on;
    T_render= toc(tm1);
    fprintf('%d-th tree rendering done, time: %f\n', i, T_render);
end

set(fig,'Position', [100 100 1000 600]);
set(gca,'position',[0.0,0.05,0.99,0.95] );
set(fig,'Color', [1, 1, 1],'InvertHardcopy','off'); 
% print(gcf, '-depsc', '-r1500', './GeoNeuro_chen_No32and41')
% % Save as .pdf
% set(fig,'Units','Inches');
% pos = get(fig,'Position');
% set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(fig, './Geod_compBotanTrees', '-dpdf', '-r1000')



% Title: Geod info
S1 = sprintf('Geodesic, BotanTree %d with %d, --- ', idx1, idx2);
Lambda_str = ['\lambda_{m}=', num2str(lam_m), ', \lambda_{s}=', num2str(lam_s), ', \lambda_{p}=', num2str(lam_p), ' --- E-total=', num2str(G.E)];
EandTime_str = ['Time :', 'T_{PadAndAlign}=', num2str(T1), 's',', T_{Geod}=', num2str(T2),'s'];


title({[S1, Lambda_str] , EandTime_str, ['time=', num2str(time)]});



% Save as .png
% set(gca,'looseInset',[0 0 0 0])
print(fig,['Geod_compBotanTrees-10timeTest-', num2str(time)],'-dpng','-r600')
% saveas(fig,'results/eg_singleShape_singleCol.png')


