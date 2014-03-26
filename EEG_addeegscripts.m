function EEG_addeegscripts()
addpath(genpath('/home/cui/scripts/EEGLABv11/eeglab11_0_2_1b'));
P=mfilename('fullpath');
P1=mfilename;
script_path=P(1:(length(P)-length(P1)));
addpath(genpath(script_path));

excluded_paths={'archive','debug'};
for n = 1:length(excluded_paths)
    rmpath(genpath([script_path,excluded_paths{n}]));
end
