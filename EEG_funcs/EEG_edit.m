function EEG_edit()
P=mfilename('fullpath');
[file_path,~,~]=fileparts(P);
addpath(file_path);
edit EEG_CreateModules.m EEG_CreateParams.m EEG_PreProc.m EEG_PostProc.m
end