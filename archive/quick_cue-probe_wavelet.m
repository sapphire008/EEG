%% Jong's new POP Data Analysis

%% Set Paths
addpath(genpath('/home/gomes/EEGLABv11/'));

%% Set Group Information
subject_group = 'pilot_';
elecs_dir = '/home/cui/eeg_4pop/test_1/elec_files/';
%preprocdir = ['/home/cui/eeg_4pop/test_1/',subject_group, '/'];
preprocdir = ['/nfs/jong_exp/EEG_4POP/subjects/'];
phase = {'Cue'};

%% Extract Conditions
master_array = []

subjectID = 'pilot_991';
name_set = [subjectID,'-Analysis-Ready.set']

%cd([preprocdir,subjectID])

for i = 1:length(phase)
EEG = pop_loadset(name_set, [preprocdir, subjectID, '/', phase{i}, '/'])

cd(phase{i})
mkdir('Conditions')
cd('Conditions')

cue_green_idx = []
cue_red_idx = []

% Extract current condition indeces JUST CUE
for ep = 1:1:length(EEG.epoch)
    %current_condition is the event code formatted as double
    if isnumeric(EEG.epoch(ep).eventtype) == 1
        current_condition = EEG.epoch(ep).eventtype;
    end
    if iscell(EEG.epoch(ep).eventtype) == 1
        current_condition = EEG.epoch(ep).eventtype{1,1} ;       
    end
    if ischar(EEG.epoch(ep).eventtype) == 1
        current_condition = str2num(EEG.epoch(ep).eventtype);
    end
    
    if current_condition == 5
        cue_green_idx(end+1) = ep %storing what epoch #s belong to this condition
    end
    if current_condition == 6
        cue_red_idx(end+1) = ep
    end
    
end

% green cue
EEG = pop_select( EEG, 'trial', cue_green_idx);
EEG = pop_saveset( EEG,[subjectID,'-cue_green_all.set'])

clear EEG

% red cue
cd([preprocdir,subjectID,'/',phase{i}])
EEG = pop_loadset(name_set)
EEG = pop_select( EEG, 'trial', cue_red_idx); 
cd('Conditions')
EEG = pop_saveset( EEG,[subjectID,'-cue_red_all.set'])

clear EEG

end


%% Time-Frequency Analysis CUE (ink-word)
%NEED TO MODIFY: use time intervals for modaf
%                need to increase time resolution
%                need to modify to loop phase

%set save name (should change for additional subjects)
ERSP_save_str = 'ERSP_cue_tony_30-80Hz.mat'
ITC_save_str = 'ITC_cue_tony_30-80Hz.mat'

phase = 'Cue'

condition_type = [{'cue_red_all'} {'cue_green_all'}]

for s = [991]
    
        
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

cd([preprocdir,subjectID,'/',phase,'/Conditions/'])

cue_red_data = [];
cue_green_data=[];

EEG = pop_loadset( [subjectID,'-',condition_type{1},'.set']);

cue_red_data = EEG.data;
clear EEG

EEG = pop_loadset( [subjectID,'-',condition_type{2},'.set']);

cue_green_data = EEG.data;

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
 
 [tmp_ersp itc powbase times freqs_new]=newtimef({cue_red_data(n,:,:) cue_green_data(n,:,:)}, ...
              EEG.pnts,[EEG.xmin EEG.xmax]*1000, EEG.srate,[6 1] ,'timesout',[-400:5:2700],'freqs',[30:1:80], ...
              'elocs', EEG.chanlocs,'baseline',[-200 0],'chaninfo', EEG.chaninfo,...
              'padratio', 16, 'plotphase','off','plotersp','off','plotitc','off','verbose','off');
 
 
        eval_lead_incong = ['ERSP.','cue','{',num2str(s),'}.incong.lead{',num2str(elec),'}=tmp_ersp{1,1};']
        eval(eval_lead_incong);

        eval_lead_cong = ['ERSP.','cue','{',num2str(s),'}.cong.lead{',num2str(elec),'}=tmp_ersp{1,2};']
        eval(eval_lead_cong);

        eval_lead_contrast = ['ERSP.','cue','{',num2str(s),'}.contrast.lead{',num2str(elec),'}=tmp_ersp{1,3};']
        eval(eval_lead_contrast);

       
        eval_lead_incong = ['ITC.','cue','{',num2str(s),'}.incong.lead{',num2str(elec),'}=itc{1,1};']
        eval(eval_lead_incong);

        eval_lead_cong = ['ITC.','cue','{',num2str(s),'}.cong.lead{',num2str(elec),'}=itc{1,2};']
        eval(eval_lead_cong);

        eval_lead_contrast = ['ITC.','cue','{',num2str(s),'}.contrast.lead{',num2str(elec),'}=itc{1,3};']
        eval(eval_lead_contrast);
                
end

end

%end

    freqs = freqs_new
    timesout=times
    mkdir([preprocdir,'/',phase,'/Data/'])
    cd([preprocdir,'/',phase,'/Data/'])
    save(ERSP_save_str, 'ERSP', 'timesout', 'freqs')
    save(ITC_save_str, 'ITC', 'timesout', 'freqs')
    
    %% Load in relavent ERSP

%if cue
condition_type = [{'cue'}]; 
trial_type = [{'incong'} {'cong'} {'contrast'}]

PFC_type = [{'all_PFC'} {'left_PFC'} {'mid_PFC'} {'right_PFC'}]
%PFC_type = [{'all_Central'} {'left_Central'} {'mid_Central'} {'right_Central'}]

%%get structs

for l = 1:length(PFC_type)
    
    for t = 1:length(trial_type)
   
      
for k= 1:length(condition_type)

    cond = condition_type{k}
    
j=1;
cond_struct_all = eval(['ERSP.',cond]);
for i=1:length(cond_struct_all)
    if length(cond_struct_all{i}) ~= 0
        eval([cond,'_struct{j} = cond_struct_all{i}'])
        j = j+1
    end
end

clear i j cond_struct_all



%yield <subject X electrode X freqency X time_interval>

for i=1:length(eval([cond,'_struct']))
    
        PFC = []
        if isequal(PFC_type{l}, 'all_PFC')
            PFC = [15 16 17 6 7 8 4 1 2 9 10 11 18 19 20 5 12 13 14 21 22 23 ]
        end
        if isequal(PFC_type{l}, 'left_PFC')
            PFC = [15 16 17 6 7 8 4 1]
        end
        if isequal(PFC_type{l}, 'mid_PFC')
            PFC = [2 9 10 11 18 19 20] 
        end
        if isequal(PFC_type{l}, 'right_PFC')
            PFC = [5 12 13 14 21 22 23] 
        end
        if isequal(PFC_type{l}, 'all_Central')
            PFC = [24 25 26 34 35 36 27 28 29 37 38 39 30 31 32 40 41 42]
        end
        if isequal(PFC_type{l}, 'left_Central')
            PFC = [24 25 26 34 35 36]
        end
        if isequal(PFC_type{l}, 'mid_Central')
            PFC = [27 28 29 37 38 39]            
        end
        if isequal(PFC_type{l}, 'right_Central')
            PFC = [30 31 32 40 41 42]
        end
            
      for elec = 1:length(PFC)
        contrast_leads = eval([cond,'_struct{i}.',trial_type{t},'.lead']);
       
       if isempty(contrast_leads{PFC(elec)}) == 0
                gamma_struct(i,elec,:,:) = contrast_leads{PFC(elec)}(:,:);
        else
                gamma_struct(i,elec,1:length(freqs),1:length(timesout)) = NaN;
        end
        
    end
          
end

%squeeze the mean across electrode_type; yields <subject X frequency X
%time_interval>. squeeze again for mean <frequency X time_interval>

%mean_struct = squeeze(nanmean(squeeze(nanmean(gamma_struct,2))));
mean_struct = (squeeze(nanmean(gamma_struct,2)));

eval([PFC_type{l},'_',cond,'_',trial_type{t},'_struct = mean_struct;']);


clearvars i k gamma_struct mean_struct
end  
 
end

clear PFC gamma_struct
end

clearvars -except *struct* freqs timesout ERSP

%% save spectrograms 

group =['sc']

%if cue
% condition_type = [{'cue'}]; 
% trial_type = 'ink_word' %one at a time


%if cue sort
condition_type = [{'cue'}]; 
trial_type = 'contrast'


PFC_type = [{'all_PFC'} {'left_PFC'} {'mid_PFC'} {'right_PFC'}]

for k = 1:length(condition_type)
    
for j = 1: length(PFC_type)
    
    sjset = eval([PFC_type{j},'_',condition_type{k},'_',trial_type,'_struct'])

    close all

plot_name = [group,'_',condition_type{k},'_',PFC_type{j}]

h = tftopo(sjset,timesout,freqs,'title',[plot_name],'verbose','off');
    caxis([-1.5 1.5])
    
    set(gca,'YTick',[freqs(1):5:freqs(end)])
    
    colorbar
    figure;
    saveas(1, [plot_name,'.tif'],'tiff')
end
end

    