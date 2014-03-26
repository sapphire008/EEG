%% Set Paths
addpath(genpath('/home/gomes/EEGLABv11/'));

%% Set Group Information
subject_group = 'npt';
phase = {'Cue','Probe'}
elecs_dir = '/nfs/erp-pop/Jong_newPOP/elec_files/';
preprocdir = ['/nfs/erp-pop/Jong_newPOP/',subject_group,'/'];
dirRawEEG = '/nfs/erp-data/Jong 4POP/';
%% Load in Subjects 
% (here, make sure they were manually loaded into EEGLAB and saved as .set

for s = [1] %sc is 16:24,26,99
    
    % Create Subject ID
    if s<100
        zero_pad = '0'
    end
    if s <100 && s<10
        zero_pad = '00'
    end
    if s>=100
        zero_pad = ''
    end
    
    subjectID = [subject_group,zero_pad,num2str(s)]
    
    EEG = pop_loadcnt([dirRawEEG, subjectID, '.cnt'], 'dataformat', 'int32', 'keystroke', 'on');

    %move into directory
    cd(preprocdir)
    mkdir(subjectID);
    cd(subjectID)
    
    %save unepocked set
    EEG = pop_saveset( EEG,[subjectID,'-unepoched.set'])
    
    % Load unepoched EEG data
    %EEG = pop_loadset( [subjectID,'-unepoched.set'])
    
    %Change cues and probes for correct events
    
    % load electrode locations
    EEG.chanlocs=readlocs([elecs_dir,'/Sixty_Four_channel_neuroscan_cap.ced']);
    %load bad electrodes
    eval(['load ',elecs_dir,'/bad_elecs_npt.mat'])
    %remove VEO HEO EKG EMG (= 65-68) CB1 CB2 (= 64,60) M1 M2 (= 33 43)
    EEG = pop_select( EEG, 'nochannel', [33 43 60 64 65 66 67 68]); 
    
    %remove listed electrodes
    bad_elecs_label = bad_elecs{s,2}
    EEG = pop_select( EEG, 'nochannel', bad_elecs_label); 
    
    % Downsample data
    EEG = pop_resample(EEG,250);

    % high-pass filter at 0.5 Hz
    EEG.data = eegfilt(EEG.data, EEG.srate, .5, 0);

    % Save Pre-epoched Data
    EEG = pop_saveset( EEG,[subjectID,'-pre-epoched.set'])
    
 % Epoch cue and probe
    for i=1:length(phase)
        clear EEG;
        
        EEG = pop_loadset( [subjectID,'-pre-epoched.set'])
        if (isequal(phase{i}, 'Cue'))
            EEG = pop_epoch( EEG, {5 6}, [-0.4 1.7], 'epochinfo', 'yes');
        elseif (isequal(phase{i}, 'Probe'))
            EEG = pop_epoch( EEG, {20 21 22 23 30 31 32 33}, [-0.4 1.2], 'epochinfo', 'yes');
        end
    
        % Use average reference
        EEG = pop_reref( EEG, [], 'refstate',0);

        % baseline subtract
        EEG = pop_rmbase( EEG, [-400 0]);

        mkdir(phase{i});

        % save the epoched data set
        EEG = pop_saveset( EEG,[subjectID,'-epoched.set'], [preprocdir, subjectID, '/', phase{i}, '/'])

        % Artifact reject

        % reject epochs based on probability
        [EEG, locthresh, globthresh, nrej] = pop_jointprob(EEG,1,[1:EEG.nbchan] ,3,5,1,0,0);

        %Reject marked epochs and save the dataset
        EEG = eeg_rejsuperpose( EEG, 1, 1, 1, 1, 1, 1, 1, 1);
        EEG = pop_rejepoch( EEG, find(EEG.reject.rejglobal), 0);

        % Reject based on voltage deflection
        %EEG = pop_eegthresh(EEG,1,[1:EEG.nbchan],-150,150,EEG.xmin,EEG.xmax,1,0);

        % save dataset
        EEG = pop_saveset( EEG,[subjectID,'-epoched-PostR.set'], [preprocdir, subjectID, '/', phase{i}, '/'])


        clear EEG locthresh globthresh nrej bad_elecs
    end
end

%% ICA
%runs both cue and probe
for s =[1]
    % Create Subject ID
    if s<100
      zero_pad = '0'
    end
    if s <100 && s<10
            zero_pad = '00'
    end
    if s>=100
        zero_pad = ''
    end
    
    subjectID = [subject_group,zero_pad,num2str(s)]
    cd([preprocdir,subjectID])

    for i = 1:length(phase)
        EEG = pop_loadset( [subjectID,'-epoched-PostR.set'], [preprocdir, subjectID, '/', phase{i}, '/'])

        %for all elecs
        EEG = pop_runica( EEG, 'icatype', 'runica','extended',1);

        EEG = pop_saveset( EEG,[subjectID,'-ICA-ready.set'], [preprocdir, subjectID, '/', phase{i}, '/'])

        clear EEG
    end
end

%% View ICs

%USER-DEF VARS
s = [1]
myPhase = 'Cue'


% Create Subject ID
if s<100
  zero_pad = '0'
end
if s <100 && s<10
        zero_pad = '00'
end
if s>=100
    zero_pad = ''
end
    
subjectID = [subject_group,zero_pad,num2str(s)]

cd([preprocdir,subjectID, '/',myPhase])

clear EEG
EEG = pop_loadset([subjectID,'-ICA-ready.set'])

close all;

% for all elecs
pop_prop(EEG,0,1:15,[])

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
for s =[1]
    % Create Subject ID
    if s<100
      zero_pad = '0'
    end
    if s <100 && s<10
            zero_pad = '00'
    end
    if s>=100
        zero_pad = ''
    end

    subjectID = [subject_group,zero_pad,num2str(s)]

    cd([preprocdir,subjectID])
    
    for i = 1:length(phase)
        EEG = pop_loadset( [subjectID,'-ICA-ready.set'], [preprocdir, subjectID, '/', phase{i}, '/'])

        eval(['load ',preprocdir,phase{i},'_removed_components.mat'])
        remove_comps = []
        remove_comps = removed_components{s,2}

        EEG = pop_subcomp(EEG,remove_comps,0);

        EEG = pop_saveset( EEG,[subjectID,'-Post-ICA.set'], [preprocdir, subjectID, '/', phase{i}, '/'])

        clear EEG
    end
end

%% Final Artifact Reject

for s =[0]
    % Create Subject ID
    if s<100
      zero_pad = '0'
    end
    if s <100 && s<10
            zero_pad = '00'
    end
    if s>=100
        zero_pad = ''
    end

    subjectID = [subject_group,zero_pad,num2str(s)]

    cd([preprocdir,subjectID])
    
    for i = 1:length(phase)
        EEG = pop_loadset( [subjectID,'-Post-ICA.set'], [preprocdir, subjectID, '/', phase{i}, '/'])

        EEG = pop_eegthresh(EEG,1,[1:EEG.nbchan],-50,50,EEG.xmin,EEG.xmax,1,0);

        EEG = pop_saveset( EEG,[subjectID,'-Analysis-Ready.set'], [preprocdir, subjectID, '/', phase{i}, '/'])
    end
end