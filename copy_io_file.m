clear all; close all; clc

prefix = 'Y:\Kiah\data\R4';
project_folder_all = dir(prefix);
source = 'Y:\Kiah\data\R4\source_io_yaml_file\io.yaml';

for k = 1:numel(project_folder_all)
    destination = [prefix,filesep,project_folder_all(k).name];
    if exist([destination,'\io.yaml'],'file')
        % make a new folder - old yaml file
        old_yaml_k = [destination,filesep,'old yaml file'];
        mkdir(old_yaml_k)
        
        % move old yaml file to that folder
        movefile([destination,'\io.yaml'],old_yaml_k);
    end
        
    if contains(destination,'ToneWater') == 1
        copyfile(source,destination);
    end 
    k
end
