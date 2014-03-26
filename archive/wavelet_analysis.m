%% Jong's new POP Data Analysis

%% Set Paths
addpath(genpath('/home/gomes/EEGLABv11/'));

%% Set Group Information
% subject_group = 'npt';
% elecs_dir = '/nfs/erp-pop/Jong_newPOP/elec_files/';
% preprocdir = ['/nfserp-pop/Jong_newPOP/',subject_group, '/'];
% phase = {'Cue','Probe'};

subject_group = 'sc';
elecs_dir = '/home/cui/eeg_4pop/test_1/elec_files/';
preprocdir = ['/nfs/erp-pop/Jong_newPOP/',subject_group, '/'];
phase = {'Cue','Probe'};


%% Extract Conditions
master_array = []

for s = [16]
    
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

    name_set = [subjectID,'-Analysis-Ready.set']

cd([preprocdir,subjectID])

for i = 1:length(phase)
EEG = pop_loadset(name_set, [preprocdir, subjectID, '/', phase{i}, '/'])

mkdir('Conditions')
cd('Conditions')
    
ink_cue_idx = []
word_cue_idx = []
word_incong_event = [10 20 30 40 27]
word_incong_probe_idx = []
word_cong_event = [11 21 31 41 17]
word_cong_probe_idx = []
ink_incong_event = [12 22 32 42 37]
ink_incong_probe_idx = []
ink_cong_event = [13 23 33 43 47]
ink_cong_probe_idx = []
current_condition = []
total_epoch = size(EEG.data,3)

% Extract current condition indeces
for ep = 1:1:length(EEG.epoch)
    %current_condition is the event code formatted as double
    if isnumeric(EEG.epoch(1).eventtype) == 1
        current_condition = EEG.epoch(ep).eventtype
    end
    if iscell(EEG.epoch(1).eventtype) == 1
        current_condition = EEG.epoch(ep).eventtype{1,1}        
    end
    if ischar(current_condition) == 1
        current_condition = str2num(current_condition)
    end
    
    if current_condition == 55
        ink_cue_idx(end+1) = ep %storing what epoch #s belong to this condition
    end
    if current_condition == 45
        word_cue_idx(end+1) = ep
    end
    
end




% word cue
EEG = pop_select( EEG, 'trial', word_cue_idx);
EEG = pop_saveset( EEG,[subjectID,'-word_cue.set'])

clear EEG

% ink cue
cd([preprocdir,subjectID])
EEG = pop_loadset(name_set)
EEG = pop_select( EEG, 'trial', ink_cue_idx); 
cd('Conditions')
EEG = pop_saveset( EEG,[subjectID,'-ink_cue.set'])

clear EEG

% word incong probe
cd([preprocdir,subjectID])
EEG = pop_loadset(name_set)
EEG = pop_select( EEG, 'trial', word_incong_probe_idx); 
cd('Conditions')
EEG = pop_saveset( EEG,[subjectID,'-word_incong_probe_all.set'])

clear EEG

% word cong probe
cd([preprocdir,subjectID])
EEG = pop_loadset(name_set)
EEG = pop_select( EEG, 'trial', word_cong_probe_idx); 
cd('Conditions')
EEG = pop_saveset( EEG,[subjectID,'-word_cong_probe_all.set'])

clear EEG

% ink incong probe
cd([preprocdir,subjectID])
EEG = pop_loadset(name_set)
EEG = pop_select( EEG, 'trial', ink_incong_probe_idx); 
cd('Conditions')
EEG = pop_saveset( EEG,[subjectID,'-ink_incong_probe_all.set'])

clear EEG

% ink cong probe
cd([preprocdir,subjectID])
EEG = pop_loadset(name_set)
EEG = pop_select( EEG, 'trial', ink_cong_probe_idx); 
cd('Conditions')
EEG = pop_saveset( EEG,[subjectID,'-ink_cong_probe_all.set'])

% all cong probe
cd([preprocdir,subjectID])
EEG = pop_loadset(name_set)
EEG = pop_select( EEG, 'trial', [ink_cong_probe_idx word_cong_probe_idx]); 
cd('Conditions')
EEG = pop_saveset( EEG,[subjectID,'-all_cong_probe.set'])

% all incong probe
cd([preprocdir,subjectID])
EEG = pop_loadset(name_set)
EEG = pop_select( EEG, 'trial', [ink_incong_probe_idx word_incong_probe_idx]); 
cd('Conditions')
EEG = pop_saveset( EEG,[subjectID,'-all_incong_probe.set'])


clear EEG

master_array(end+1,1) = s % subject
master_array(end,2) = total_epoch % total epochs
master_array(end,3) = length(word_cue_idx) % number of word cue
master_array(end,4) = length(ink_cue_idx) % number of ink cue
master_array(end,5) = length(word_incong_probe_idx) % number of word incong probe
master_array(end,6) = length(word_cong_probe_idx) % number of word cong probe
master_array(end,7) = length(ink_incong_probe_idx) % number of ink incong probe
master_array(end,8) = length(ink_cong_probe_idx) % number of ink cong probe

end


cd([preprocdir,'Data'])
xlswrite([subject_group,'_trial-count.xls'], master_array)
end
%% Time-Frequency Analysis CUE (ink-word)
%NEED TO MODIFY: use time intervals for modaf
%                need to increase time resolution
%                need to modify to loop phase

%set save name (should change for additional subjects)
ERSP_save_str = 'ERSP_cue_36.mat'
ITC_save_str = 'ITC_cue_36.mat'

condition_type = [{'ink_cue'} {'word_cue'}]

for s = [36]
    
        
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

cd([preprocdir,subjectID,'/Conditions/'])

word_cue_data = [];
ink_cue_data=[];

EEG = pop_loadset( [subjectID,'-',condition_type{2},'.set']);

word_cue_data = EEG.data;
clear EEG

EEG = pop_loadset( [subjectID,'-',condition_type{1},'.set']);

ink_cue_data = EEG.data;

elecs = {EEG.chanlocs.type}

  
for n = 1:length(elecs)
        elec = elecs{n}

 if ischar(elec) == 1
        elec = str2num(elec);
 end
 
 tmp_ersp = []
 itc = []
 
 clc
 disp(['Subject: ', num2str(s), ' | Electrode: ', num2str(n),'/',num2str(length(elecs))])
 
 [tmp_ersp itc powbase times freqs_new]=newtimef({ink_cue_data(n,:,:) word_cue_data(n,:,:)}, ...
              EEG.pnts,[EEG.xmin EEG.xmax]*1000, EEG.srate,[2 6] ,'timesout',[-1000:4:3500],'freqs',[2:1:80], ...
              'elocs', EEG.chanlocs,'baseline',[-400 0],'chaninfo', EEG.chaninfo,...
              'padratio', 16, 'plotphase','off','plotersp','off','plotitc','off','verbose','off');
 
 
        eval_lead_incong = ['ERSP.','cue','{',num2str(s),'}.ink.lead{',num2str(elec),'}=tmp_ersp{1,1};']
        eval(eval_lead_incong);

        eval_lead_cong = ['ERSP.','cue','{',num2str(s),'}.word.lead{',num2str(elec),'}=tmp_ersp{1,2};']
        eval(eval_lead_cong);

        eval_lead_contrast = ['ERSP.','cue','{',num2str(s),'}.ink_word.lead{',num2str(elec),'}=tmp_ersp{1,3};']
        eval(eval_lead_contrast);

       
        eval_lead_incong = ['ITC.','cue','{',num2str(s),'}.ink.lead{',num2str(elec),'}=itc{1,1};']
        eval(eval_lead_incong);

        eval_lead_cong = ['ITC.','cue','{',num2str(s),'}.word.lead{',num2str(elec),'}=itc{1,2};']
        eval(eval_lead_cong);

        eval_lead_contrast = ['ITC.','cue','{',num2str(s),'}.ink_word.lead{',num2str(elec),'}=itc{1,3};']
        eval(eval_lead_contrast);
                
end

end

%end

    freqs = freqs_new
    timesout=times
    cd([preprocdir,'/Data/'])
    save(ERSP_save_str, 'ERSP', 'timesout', 'freqs')
    save(ITC_save_str, 'ITC', 'timesout', 'freqs')
    