function EEG_PreProc(work_dir)
%make sure current working directory contains workspace.mat
if nargin<1
    work_dir=pwd;
elseif strcmpi(work_dir(end),'/')
    work_dir=work_dir(1:end-1);
end
load([work_dir,'/workspace.mat']);
%Add EEG processing script to the search directory
addpath(genpath(params.Script_dir));
%Add EEGLAB package to the search directory
EEG_addeegscripts;
%diary([params.subjectID '_log.txt']);%record logs during data processing
%% Start Preprocessing
%a). Convert from .cnt to .set format (-unepoched.set)
EEG = pop_loadcnt([params.RawEEG_dir,'/',...
    getfield(dir([params.RawEEG_dir,'/*.cnt']),'name')],...%get filename
    'dataformat', 'int32', 'keystroke', 'on');%load raw cnt file
%save unepoched set
EEG = pop_saveset(EEG,'filepath',params.Save_dir,...
    'filename',[params.subjectID,'-unepoched.set']);

%b). Remove bad electrodes according to documentation (-pre-epoched.set)
% load electrode locations
EEG.chanlocs=readlocs([elecs.directory,elecs.map]); %#ok<NODEF>
%remove default and documented bad channels
EEG = pop_select(EEG,'nochannel',elecs.default_bad);
EEG = pop_select(EEG,'nochannel',elecs.listBad);

%c). Downsample data to 250 Hz
EEG = pop_resample(EEG,params.downSampleRate);

%d). High-pass filter at 0.5 Hz
EEG.data = eegfilt(EEG.data, EEG.srate, params.highPassFilt, 0,0,[],0,...
    'fir1',0);

%e). Save Pre-epoched Data
EEG = pop_saveset(EEG,'filepath',params.Save_dir,...
    'filename',[params.subjectID,'-pre-epoched.set']);

%% f). Epoch according to behave condition (-epoched.set)
for n=1:numPhase%for each phase of the task
    disp(['%%%%%%%%%%%% Current Phase: ',Phase{n,1},' %%%%%%%%%%%%']);
    
    clear EEG locthresh globthresh nrej bad_elecs;
    EEG = pop_loadset('filepath',[params.Save_dir ...
        '/'],'filename',[params.subjectID,'-pre-epoched.set']);
    EEG=pop_epoch(EEG,Phase{n,2},Phase{n,3}(1,[1,end])+Phase{n,4},...
        'epochinfo','yes');

    % Use average reference
    EEG = pop_reref(EEG, [], 'refstate',0);
    % baseline subtract
    switch upper(Phase{n,5})
        case {'ON'}%if marker inserted at onset
            EEG = pop_rmbase(EEG, [-params.baseline_leng 0]);
        case {'OFF'}%if marker inserted at offset
            disp('Removing baseline by conditions...');
            EEG = EEG_rmbase(EEG,Task,Condition,Phase(n,:),...
                params.baseline_leng);
    end
    %Check for remaining broken electrodes by inspecting time series
    [EEG,elecs.remBad.(Phase{n,1})]=EEG_Rejector(EEG,1);
    
     %make a directory for the current phase
    mkdir(Phase{n,1});
     % save the epoched data set
    EEG = pop_saveset(EEG,'filepath',[params.Save_dir, Phase{n,1}, ...
        '/'],'filename',[params.subjectID,'-epoched.set']);
    save([params.Save_dir,'workspace.mat'],'elecs','-append');
end

%% h). Artifact Reject (-epoched-PostR.set) and ICA (-ICA-ready.set)
clear n;
for n = 1:numPhase
    EEG = pop_loadset('filepath',[params.Save_dir Phase{n,1} ...
        '/'],'filename',[params.subjectID,'-epoched.set']);
    % Artifact reject
    % case 1: reject epochs based on probability
    EEG=EEG_artifactRejection(EEG,1);
    % case 2: reject epochs based voltage deflection
    %leaving only voltage between [-150 150] mV
    %EEG=EEG_artifactRejection(EEG,2,150);
   
    EEG = pop_saveset(EEG,...
        'filepath',[params.Save_dir Phase{n,1} '/'],...
        'filename',[params.subjectID,'-epoched-PostR.set']); 
    
    %run ICA for all elecs
    EEG = pop_runica(EEG, 'icatype', 'runica','extended',1);
    %save ICA result
    EEG = pop_saveset(EEG,...
        'filepath',[params.Save_dir Phase{n,1} '/'],...
        'filename',[params.subjectID,'-ICA-ready.set']);

    clear EEG;
end

%% i) View&Reject ICs (-Post-ICA),Final Artifact Reject (-Analysis-Ready)
clear n;
for n=1:numPhase
    EEG = pop_loadset('filepath',[params.Save_dir Phase{n,1} ...
        '/'],'filename',[params.subjectID,'-ICA-ready.set']);
    %View ICA and save after rejecting components (-Post-ICA.set)
    clc;%clear command window
    diary([params.Save_dir,Phase{n,1},'/',Phase{n,1},'_ICA_',date,'.txt']);
    disp(['%%%%%%%%%%%% Current Phase: ',Phase{n,1},' %%%%%%%%%%%%']);
    %reject the ICA components
    [EEG,bad_comps.(Phase{n,1})]=EEG_Rejector(EEG,2,params,'numComp',20,...
        'PlotOff',0);
    %save to complete IC rejection
    EEG = pop_saveset(EEG,...
        'filepath',[params.Save_dir Phase{n,1} '/'],...
        'filename',[params.subjectID,'-Post-ICA.set']);

    % Final artifact rejection (-Analysis-Ready.set)
    % case 2: reject epochs based voltage deflection
    EEG=EEG_artifactRejection(EEG,2,'thresh',params.artThresh);
    % Double check if the epoches are labeled with correct markers
    EEG=EEG_checkEpoch(EEG,'InvalidMarker',255,'DisplayEpoch',0);
    %save after final artifact reject and marker check
    EEG = pop_saveset(EEG,'filepath',[params.Save_dir Phase{n,1} ...
        '/'],'filename',[params.subjectID,'-Analysis-Ready.set']);
    numEpoch.(Phase{n,1})=length(EEG.epoch);
    diary off;
    save([params.Save_dir,'workspace.mat'],'elecs','bad_comps',...
    'numEpoch','-append');
end

end