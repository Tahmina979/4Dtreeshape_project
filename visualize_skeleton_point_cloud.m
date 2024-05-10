clc
clear
clear all;

addpath('dataset')
% #### DIRECTORY OF THE DATASET #####

path_to_swc="address of the directory of .swc files in dataset folder"
path=dir(strcat(path_to_swc,"*.swc"));

for i=1:length(path)
    file = fopen(strcat(path_to_swc,path(i).name));
    data = zeros(0,3);
    while(~feof(file))
    fline = fgetl(file);  
    if fline(1) == '#'
        continue;
    end    
    A = sscanf(fline, '%f %f %f %f %f %f %f');
    A1=A(3:5);
    data = [data; A1']
    end
    fclose(file);
    newObj=pointCloud(data); 
    figure(1);
    subplot(1, length(path), i);
    pcshow(newObj);
end
return;
