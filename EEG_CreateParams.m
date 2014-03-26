function EEG_CreateParams(subjID,doc_bad_elecs,task_Name)
%EEG_CreateParams(subjectID,documented_bad_elecs,task_Name)
%Create parameters necessary for processing
%subjectID: a string
%documented_bad_elecs: cell array of strings
%make sure the folder name and the file names under the folder use
%subjectID
%task_Name: name of the task currently processing
%   select among: 
%           '4POP'
%           'WM'
% Optional input:
%   'PhaseType' : 'simple'  --load simple phase type
%                 'complex' --load complex phase type (default)
%% open the .m script if there is no input
if isempty(subjID) || isempty(doc_bad_elecs)
    P=mfilename('fullpath');
    edit(P);
end

%% Set subject related parameters   
%Subject Profile
%subject ID: Make sure the file name and folder name matches this
params.subjectID=subjID;
%documented bad electrodes for this subject
elecs.listBad=doc_bad_elecs;
%base directory of the folder structures
base_dir='/nfs/jong_exp/EEG_PFC/';
%% ---The following  do not always change from subject to subject---
%% a) Other subject profiles
%raw data directory: containing .cnt EEG raw data, and .edat & .txt
%behavior data

params.RawEEG_dir=[base_dir,'subjects/',params.subjectID,'/',...
    task_Name,'/'];
%save directory: also saves the parameters as a .mat file to used later
params.Save_dir=[base_dir,'subjects/',params.subjectID,'/',task_Name,'/'];
%script directory
params.Script_dir=[base_dir,'scripts/'];
%module directory
params.Module_dir=[base_dir, 'scripts/modules/'];
%.txt-->.csv Pearl script directory
%txt2csv_script_dir='/nfs/jong_exp/EEG_PFC/scripts/edat2csv_SC_v1.pl';

%% b) Setting the current task conditions
%load module
if ~exist([params.Module_dir,'defaultParams.mat'],'file')
    flag=input('defaultParams.mat does not exist, create one? (Y/N):','s');
    flag=upper(flag);
    switch flag
        case {'Y'}
            EEG_CreateModules();
            load([params.Module_dir,'defaultParams.mat']);
        case {'N'}
            edit('EEG_CreateModules.m');
            return;
    end
else
    load([params.Module_dir,'defaultParams.mat']);
end

%make task structure from default parameters
if strcmpi(task_Name,'4POP') %incase of 4POP, which starts with a number
    task_Name='FourPOP';
end
Task=taskSet.(task_Name).Name;%task name
Phase=taskSet.(task_Name).Phase;%set the behavioral task: 4POP
%for the selected task, calculating behavioral task profile
numPhase=length(Phase(:,1));%number of conditions
Condition=taskSet.(task_Name).Conditions;%set the conditions in the task
%numCond=length(Condition(:,1));%calculate number of conditions
Accuracy=taskSet.(task_Name).Accuracy;

%% c). Electrode profile
elecs.directory = [base_dir,'scripts/modules/'];
elecs.map='Sixty_Four_channel_neuroscan_cap.ced';
%remove VEO HEO EKG EMG (= 65-68) CB1 CB2 (= 64,60) M1 M2 (= 33 43)
elecs.default_bad=[33 43 60 64 65 66 67 68];%default bad channel
elecs.numChan=64;%number of channels of the cap used

%% d). Misc. parameters used during processing
params.downSampleRate=250;%down sample rate[Hz]
params.highPassFilt=0.5;%high pass filter frequency[Hz]
params.baseline_leng=400;%time before a phase used as baseline[ms]
params.artThresh=50;%threshold for artifact rejection based on voltage[mV]
params.bandInterest='gamma';%frequency band of interest

%% Save parameters for later use during processing
% behavior if workspace already existed
if exist([params.RawEEG_dir,'workspace.mat'],'file')
    inquiry=input([subjID, ...
        '''s workspace.mat already existed. Replace? (Y/N):'],'s');
    switch inquiry
        case {'Y','y'}
            eval(['!rm ',params.RawEEG_dir,'workspace.mat']);
        case {'N','n'}
            return;
    end
end
save([params.Save_dir,'workspace.mat'],'params','elecs','Task','Phase',...
    'numPhase','Condition','Accuracy');
disp('A workspace is created at:');
disp(params.Save_dir);
end