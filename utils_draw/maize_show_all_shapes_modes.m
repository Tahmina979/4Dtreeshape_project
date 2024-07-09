function maize_show_all_shapes_modes(CT,count)

addpath Yanirk_code


default_c = get(gca,'colororder');
% c = get_colormap(100);
% 
% 

for i=1:numel(CT)
   % tree=load(strcat("C:\Users\34719017\Documents\MATLAB\complexTrees_code-main\utils_data\NeuroData\sequence_1\", mat_files(i).name));
    %CT{i} =tree.compTrees;
    fprintf('Drawing %d-th/%d tree...', i, length(CT));
    tm = tic;
    bol=isempty(CT{i});
    if bol==1
       if i==6
       % text(X0(1)-18, Y0(1), Z0(1)-45, label(count), 'FontSize', 16);
       end
        continue;
    end

    if i==1
       
        x_added = 0;
  
    
     elseif i==5
        x_added = (mod(i-1, 10)) * 78;
     elseif i==6
        x_added = (mod(i-1, 10)) * 82;
     elseif i==7
        %x_added = (mod(i-1, 10)) * 106;
        x_added = (mod(i-1, 10)) * 94;
     elseif i==8
        %x_added = (mod(i-1, 10)) * 122;
        x_added = (mod(i-1, 10)) * 104;
     elseif i==9
          %x_added = (mod(i-1, 10)) * 138;
          x_added = (mod(i-1, 10)) * 116;
     elseif i==10
          x_added = (mod(i-1, 10)) * 125;
          %x_added = (mod(i-1, 10)) * 150;  activate for whole tree
     else
          x_added = (mod(i-1, 10)) * 75;
     end

      y_added = ceil(i/10) * (-80);

 
   if count==1
      
        z_added=0;
   end
   if count==2
     
      %  z_added = ceil(i/10) * (-180); % active for whole tree
        z_added = ceil(i/10) * (-400);
       
        
   end
   if count==3
      
        z_added = ceil(i/10) * (-800);
       
   end

   if count==4
      
        z_added = ceil(i/10) * (-1200);
        
        
   end

   if count==5
     
        z_added = ceil(i/10) * (-1600);
       
   
   end
   if count==6
      
        z_added = ceil(i/10) * (-2000);
       
   end
    if count==7
     
        z_added = ceil(i/10) * (-2400);
    end
    

   % CT{i}
    clear first_pt_x first_pt_y first_pt_z
    
    first_pt_x = (CT{i}.beta0(1, 1));
    first_pt_y = (CT{i}.beta0(2, 1));
    first_pt_z = (CT{i}.beta0(3, 1));
    
    X0 = CT{i}.beta0(1, :)+ x_added-first_pt_x ;
    Y0 = CT{i}.beta0(2, :) + y_added-first_pt_y ;
    Z0 = CT{i}.beta0(3, :)+z_added- first_pt_z;
   % R0 = CT{i}.beta0_rad;
   plot3(X0, Y0, Z0, 'Color',default_c(2,:), 'LineWidth', 1); hold on;  
 
   
    for j=1: numel(CT{i}.beta)
        clear X Y Z
        
    if size(unique(CT{i}.beta{j}),1)==3
        
        continue;
    end
   
       
        X = CT{i}.beta{j}(1, :)+ x_added- first_pt_x ;
        Y = CT{i}.beta{j}(2, :) + y_added- first_pt_y ;
        Z = CT{i}.beta{j}(3, :)+z_added- first_pt_z;
        %R= CT{i}.beta_rad{j};
        
        plot3(X, Y, Z, 'Color',default_c(1,:), 'LineWidth', 1); hold on;  
      %text(X(length(X)), Y(length(Y)), Z(length(Z)), num2str(j), 'FontSize', 10);
        for k = 1: numel(CT{i}.beta_children{j}.beta)
            
           clear X Y Z
            X = CT{i}.beta_children{j}.beta{k}(1, :)+ x_added- first_pt_x ;
            Y = CT{i}.beta_children{j}.beta{k}(2, :)+ y_added- first_pt_y ;
            Z = CT{i}.beta_children{j}.beta{k}(3, :)+z_added- first_pt_z;
           % R= CT{i}.beta_children{j}.beta_rad{k};
            
            plot3(X, Y, Z, 'Color',default_c(3,:), 'LineWidth', 1); hold on;
           % text(X(length(X)), Y(length(Y)), Z(length(Z)), num2str(k), 'FontSize', 10);
       
            
            for t = 1: numel(CT{i}.beta_children{j}.beta_children{k}.beta)
                clear X Y Z
                if size(unique(CT{i}.beta_children{j}.beta_children{k}.beta{t}),1)<=9
                    continue;
                end
                X = CT{i}.beta_children{j}.beta_children{k}.beta{t}(1, :)+ x_added- first_pt_x ;
                Y = CT{i}.beta_children{j}.beta_children{k}.beta{t}(2, :) + y_added- first_pt_y ;
                Z = CT{i}.beta_children{j}.beta_children{k}.beta{t}(3, :)+z_added- first_pt_z;
               % R = CT{i}.beta_children{j}.beta_children{k}.beta_rad{t};
                plot3(X, Y, Z, 'Color',default_c(4,:), 'LineWidth', 1); hold on;
             % text(X(length(X)), Y(length(Y)), Z(length(Z)), num2str(t), 'FontSize', 10);
            end
    
            
        end
        
    end
    
    T_cost = toc(tm);
    fprintf('-Done, time cost: %.2f\n', T_cost)
        
end

end