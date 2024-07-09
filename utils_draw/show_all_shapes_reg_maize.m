function show_all_shapes_reg_maize(CT,count)

addpath Yanirk_code

% c = distinguishable_colors(70);
label=["Source","Target with missing samples","Target after interpolation","Target after temporal registration","Groud truth target"];
textt=["03-13","03-15","03-17","03-19","03-21","03-24","03-25"];
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
       if i==4
           text(X0(1)-18, Y0(1), Z0(1)-45, label(count), 'FontSize', 12);
       end
       
        continue;
    end
   % y_added = ceil(i/10) * (-80);
   % x_added = (mod(i-1, 10)) * 90;
     if i==1
        x_added = 0;     
     elseif i==3
        x_added = (mod(i-1, 10)) * 100;
     elseif i==4
        x_added = (mod(i-1, 10)) * 120;
     elseif i==5
        x_added = (mod(i-1, 10)) * 150;
     elseif i==6
        x_added = (mod(i-1, 10)) * 180;
     elseif i==7
        x_added = (mod(i-1, 10)) * 210;
     else
          x_added = (mod(i-1, 10)) * 90;
     end
      y_added = ceil(i/10) * (-80);
   if count==1
        z_added=0;
   end
   if count==2  
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
        z_added = ceil(i/10) * (-900);
   end
   if count==7
        z_added = ceil(i/10) * (-630);
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
    plot3(X0, Y0, Z0, 'Color',default_c(2,:), 'LineWidth', 0.5);

   if i==4
        text(X0(1), Y0(1), Z0(1)-30, label(count), 'FontSize', 12);
   end
        text(X0(1), Y0(1)-0.1, Z0(1), textt(i), 'FontSize', 10);
    
    
    for j=1: numel(CT{i}.beta)
        clear X Y Z
        
        X = CT{i}.beta{j}(1, :)+ x_added- first_pt_x ;
        Y = CT{i}.beta{j}(2, :) + y_added- first_pt_y ;
        Z = CT{i}.beta{j}(3, :)+z_added- first_pt_z;
        %R= CT{i}.beta_rad{j};
        
        plot3(X, Y, Z, 'Color',default_c(1,:), 'LineWidth', 0.5); hold on;  
        
        for k = 1: numel(CT{i}.beta_children{j}.beta)
            clear X Y Z
            X = CT{i}.beta_children{j}.beta{k}(1, :)+ x_added- first_pt_x ;
            Y = CT{i}.beta_children{j}.beta{k}(2, :)+ y_added- first_pt_y ;
            Z = CT{i}.beta_children{j}.beta{k}(3, :)+z_added- first_pt_z;
           % R= CT{i}.beta_children{j}.beta_rad{k};
            
            plot3(X, Y, Z, 'Color',default_c(3,:), 'LineWidth', 0.5); hold on;
           % text(X(length(X)), Y(length(Y)), Z(length(Z)), num2str(k), 'FontSize', 10);
       
            
            for t = 1: numel(CT{i}.beta_children{j}.beta_children{k}.beta)
                clear X Y Z
                X = CT{i}.beta_children{j}.beta_children{k}.beta{t}(1, :)+ x_added- first_pt_x ;
                Y = CT{i}.beta_children{j}.beta_children{k}.beta{t}(2, :) + y_added- first_pt_y ;
                Z = CT{i}.beta_children{j}.beta_children{k}.beta{t}(3, :)+z_added- first_pt_z;
               % R = CT{i}.beta_children{j}.beta_children{k}.beta_rad{t};
                plot3(X, Y, Z, 'Color',default_c(4,:), 'LineWidth', 0.5); hold on;
             % text(X(length(X)), Y(length(Y)), Z(length(Z)), num2str(t), 'FontSize', 10);
            end
    
            
        end
        
    end
    
    T_cost = toc(tm);
    fprintf('-Done, time cost: %.2f\n', T_cost)
        
end

end