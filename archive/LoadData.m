%% Set Paths
addpath(genpath('/home/cui/scripts/EEGLABv11/'));

%% Set Group Information
subject_group = 'sc';
subjNum=16;%subject number 
phase = {'Cue',{5 6}, [-0.4 1.7];'Probe',{20 21 22 23 30 31 32 33},...
    [-0.4 1.2]};%conditions of the task
numPhase=length(phase(:,1));%number of conditions
elecs_dir = '/home/cui/eeg_4pop/test_1/elec_files/';
badElec_mat='bad_elecs_npt.mat';
elecs_map='Sixty_Four_channel_neuroscan_cap.ced';
preprocdir = ['/home/cui/eeg_4pop/test_1/',subject_group,'/'];
dirRawEEG='/home/cui/EEG_practice/Cued_Stroop/raw_data/';
numDig=3;%number of digits in subject number
zero_pad=repmat('0',1,(numDig-legnth(num2str(subNum))));%zero pads
subjectID = [subject_group,zero_pad,num2str(subNum)];
%remove VEO HEO EKG EMG (= 65-68) CB1 CB2 (= 64,60) M1 M2 (= 33 43)
default_bad_chan=[33 43 60 64 65 66 67 68];%default bad channel
%% Load in Subjects 
% (here, make sure they were manually loaded into EEGLAB and saved as .set

EEG = pop_loadcnt([dirRawEEG, subjectID, '-raw.cnt'], 'dataformat', ...
    'int32', 'keystroke', 'on');%load raw cnt file

%make and move to a new directory
mkdir(preprocdir,subjectID);
cd(subjectID);

%save unepoched set
EEG = pop_saveset(EEG,[subjectID,'-unepoched.set']);

% Load unepoched EEG data
%EEG = pop_loadset( [subjectID,'-unepoched.set'])

%Change cues and probes for correct events

% load electrode locations
EEG.chanlocs=readlocs([elecs_dir,'/',elecs_map]);
%load bad electrodes
load([elecs_dir,'/',badElec_mat]);

EEG = pop_select( EEG, 'nochannel', default_bad_chan); 

%remove listed electrodes
bad_elecs_label = bad_elecs{s,2};
EEG = pop_select( EEG, 'nochannel', bad_elecs_label); 

% Downsample data
EEG = pop_resample(EEG,250);

% high-pass filter at 0.5 Hz
EEG.data = eegfilt(EEG.data, EEG.srate, .5, 0);

% Save Pre-epoched Data
EEG = pop_saveset(EEG,[subjectID,'-pre-epoched.set']);

% Epoch cue and probe
clear EEG;
EEG = pop_loadset([subjectID,'-pre-epoched.set']);

for i=1:numPhase
    EEG = pop_epoch(EEG, phase{i,2}, phase{i,3}, 'epochinfo', 'yes');

    % Use average reference
    EEG = pop_reref(EEG, [], 'refstate',0);%default value 0

    % baseline subtract
    EEG = pop_rmbase(EEG, [-400 0]);%default value -400

    mkdir(phase{i,1});

    % save the epoched data set
    EEG = pop_saveset(EEG,[subjectID,'-epoched.set'], ...
        [preprocdir, subjectID, '/', phase{i,1}, '/']);

    % Artifact reject

    % reject epochs based on probability
    [EEG, locthresh, globthresh, nrej] = pop_jointprob(EEG,1,...
    [1:EEG.nbchan] ,3,5,1,0,0);

    %Reject marked epochs and save the dataset
    EEG = eeg_rejsuperpose(EEG, 1, 1, 1, 1, 1, 1, 1, 1);
    EEG = pop_rejepoch(EEG, find(EEG.reject.rejglobal), 0);

    % Reject based on voltage deflection
    %EEG = pop_eegthresh(EEG,1,[1:EEG.nbchan],-150,150,EEG.xmin,EEG.xmax,1,0);

    % save dataset
    EEG = pop_saveset( EEG,[subjectID,'-epoched-PostR.set'], ...
        [preprocdir, subjectID, '/', phase{i,1}, '/']);


    clear EEG locthresh globthresh nrej bad_elecs
end


%% ICA
%runs both cue and probe
for i = 1:numPhase
    EEG = pop_loadset( [subjectID,'-epoched-PostR.set'], ...
        [preprocdir, subjectID, '/', phase{i,1}, '/']);
    %for all elecs
    EEG = pop_runica( EEG, 'icatype', 'runica','extended',1);
    %save ICA result
    EEG = pop_saveset( EEG,[subjectID,'-ICA-ready.set'], ...
        [preprocdir, subjectID, '/', phase{i,1}, '/'])

    clear EEG
end

%% View ICs
flag='0';%default flag
while flag ~='q'
    %prompting user
    disp('Which condition would you like to see?');
    disp([num2str([1:numPhase]'),repmat('.',numPhase,1),char(phase(:,1))]);
    flag=input('Press any "q" to quit\n','s');
    try
        myPhase=phase(num2str(flag),1);%get the corresponding phase
    catch error_input
        error('Please select phase or press "q" to quit');
    end
    cd([preprocdir,subjectID, '/',myPhase,'/'])

    clear EEG
    EEG = pop_loadset([subjectID,'-ICA-ready.set'])

    close all;

    % for all elecs
    pop_prop(EEG,0,1:15,[])
end
% individual components
%figure
%eegplot(EEG.icaact,'winlength',[1], 'dispchans', [1:15]);

% 3-D models
% cd([preprocdir,num2str(group),num2str(zeroes),num2str(s)])
% EEGOUT = pop_headplot(EEG, 0, [1:15], 'Components', [3 5],'load','headplot.spl','view',[180 0],'maplimits','absmax')

%get variances
% figure
% [cmpvarorder,cmpvars,cmpframes,cmptimes,cmpsplotted,sortvar] = envtopo(EEG.data,EEG.icaweights*EEG.icasphere,'sortvar','pv','dispmaps','off')
%% Remove IC's (have "removed_components.mat ready)

cd([preprocdir,subjectID])

for i = 1:numPhase
    EEG = pop_loadset( [subjectID,'-ICA-ready.set'], [preprocdir, subjectID, '/', phase{i,1}, '/'])

    eval(['load ',preprocdir,phase{i,1},'_removed_components.mat'])
    remove_comps = []
    remove_comps = removed_components{s,2}

    EEG = pop_subcomp(EEG,remove_comps,0);

    EEG = pop_saveset( EEG,[subjectID,'-Post-ICA.set'], [preprocdir, subjectID, '/', phase{i,1}, '/'])

    clear EEG
end
%% Final Artifact Reject
cd([preprocdir,subjectID])

for i = 1:numPhase
    EEG = pop_loadset( [subjectID,'-Post-ICA.set'], [preprocdir, subjectID, '/', phase{i,1}, '/'])

    EEG = pop_eegthresh(EEG,1,[1:EEG.nbchan],-50,50,EEG.xmin,EEG.xmax,1,0);

    EEG = pop_saveset( EEG,[subjectID,'-Analysis-Ready.set'], [preprocdir, subjectID, '/', phase{i,1}, '/'])
end
