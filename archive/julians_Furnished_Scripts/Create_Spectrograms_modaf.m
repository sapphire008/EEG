%% Create Spectrograms
% Need to load ERSP into workspace first!

%% 1) load regions 

addpath(genpath('/nfs/erp-modaf/elec_files/'))
load wholehead_elecs.mat

%all PFC
all_PFC = cat(2, wholehead_elecs{1,2},wholehead_elecs{2,2}, wholehead_elecs{3,2});
   
%PFC left
left_PFC = wholehead_elecs{2,2};

%PFC mid
mid_PFC = wholehead_elecs{1,2};

%PFC right
right_PFC = wholehead_elecs{3,2};

%all Central
all_Central = cat(2, wholehead_elecs{8,2},wholehead_elecs{4,2}, wholehead_elecs{9,2});

%Central left
left_Central = wholehead_elecs{8,2};

%Central mid
mid_Central = wholehead_elecs{4,2};

%Central right
right_Central = wholehead_elecs{9,2};

%all Posterior
all_Post = cat(2, wholehead_elecs{6,2},wholehead_elecs{5,2}, wholehead_elecs{7,2});

%Posterior left
left_Post = wholehead_elecs{6,2};

%Posterior mid
mid_Post = wholehead_elecs{5,2};

%Posterior right
right_Post = wholehead_elecs{7,2};


%% 2) Form <freqsXtimes> struct
% set subject group, trial_type, and "PFC_type" (region)
% time_interval and freq_interval need to be set manually if timesout and
    % freqs don't exist with ERSP
%MAKE SURE ERSP/ITC DON'T CO-EXIST IN WORKSPACE

%NOTE: Qmoney is hard-coded as Drug contrast - Plac contrast; all other dose_types will not generate a Qmoney
    
%user-defined vars

%dose_type = [{'Plac_dose3'} {'Plac_single_d3Plac'}];
dose_type = [{'Drug'} {'Plac'}];

PFC_type = [{'all_PFC'} {'left_PFC'} {'mid_PFC'} {'right_PFC'} ...
            {'all_Central'} {'left_Central'} {'mid_Central'} {'right_Central'} ...
            {'all_Post'} {'left_Post'} {'mid_Post'} {'right_Post'}];
% PFC_type = [{'all_PFC'} {'left_PFC'} {'mid_PFC'} {'right_PFC'}];

lstConditions = [{'Green'} {'Red'} {'Contrast'}];   %Green, Red, Contrast
%lstConditions = [{'Red'}];   %Green, Red, Contrast

%subject exclusion
boolUseSubjectExclusion = false;
lstExcludedSubjects = [103,108,110,111,115];   %will not process these subjects, use to limit the MC/PC dataset to MC or PC only

%electrode-by-electrode (1/25/2013)
%This is used to diagnose if single electrodes are contributing to
%baseline activation
%Note: this will create a LOT of graphs (at most 128 per condition), so
%recommend not running too many
boolUseSingleElectrodes = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%begin algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%automatic ERSP/ITC switching
if exist('ERSP','var')
    strDataType = 'ERSP';
elseif exist('ITC','var')
    strDataType = 'ITC';
else
    error('Could not find either ERSP or ITC in workspace!');
end

time_interval = [1:length(timesout)];
freqs_interval = [1:length(freqs)];

%%get structs

for iterator = 1:length(lstConditions)
    condition = lstConditions{iterator};

    for l = 1:length(PFC_type)

        eval(['PFC =', PFC_type{l},';'])

    for k= 1:length(dose_type)

        dose = dose_type{k};
        
    

    %delete empty entries from input ERSP
    
    j=1;
    dose_struct_all = eval([strDataType,'.',dose]);
    if boolUseSubjectExclusion  %check if subject exclusion is used
        dose_struct_all(lstExcludedSubjects) = [];
        fprintf('Warning: subject exclusion is ON, excluded subjects:\t%s\n',num2str(lstExcludedSubjects))
    end
    for i=1:length(dose_struct_all)
        if ~isempty(dose_struct_all{i})
            eval([dose,'_struct{j} = dose_struct_all{i};'])
            j = j+1;
        end
    end

    clear i j dose_struct_all



    %yield <subject X electrode X freqency X time_interval>

    for i=1:length(eval([dose,'_struct']))
        for elec = 1:length(PFC)
            if isequal(condition, 'Contrast')
              contrast_leads = eval([dose,'_struct{i}.contrast.lead']);
            end
            if isequal(condition, 'Green')
                contrast_leads = eval([dose,'_struct{i}.cong.lead']);
            end
            if isequal(condition, 'Red')
                contrast_leads = eval([dose,'_struct{i}.incong.lead']);
            end
            if PFC(elec) <= size(contrast_leads,2) && isempty(contrast_leads{PFC(elec)}) == 0
                if exist('ITC','var') && ~isequal(condition, 'Contrast')    %contrast is already real numbers
                    gamma_struct(i,elec,:,:) = abs(contrast_leads{PFC(elec)}(freqs_interval,time_interval));
                else
                    gamma_struct(i,elec,:,:) = contrast_leads{PFC(elec)}(freqs_interval,time_interval);
                end
            else
                    gamma_struct(i,elec,1:length(freqs_interval),1:length(time_interval)) = NaN;
            end
        end
    end

    if ~boolUseSingleElectrodes
        %squeeze the mean across electrode_type; yields <subject X frequency X
        %time_interval>. squeeze again for mean <frequency X time_interval>

        eval([PFC_type{l},'_',dose,'_gamma_struct_initial = gamma_struct;']);

        mean_gamma_struct = squeeze(nanmean(squeeze(nanmean(gamma_struct,2))));

        eval([PFC_type{l},'_',dose,'_',condition,'_gamma_struct = mean_gamma_struct;']);
    else
        %squeeze only across subjects, keep electrode dimension

        eval([PFC_type{l},'_',dose,'_gamma_struct_initial = gamma_struct;']);

        mean_gamma_struct = squeeze(nanmean(gamma_struct,1));

        eval([PFC_type{l},'_',dose,'_',condition,'_gamma_struct = mean_gamma_struct;']);
    end

    clearvars i k gamma_struct mean_gamma_struct '*_initial'


    end

    % Take Drug minus Plac and create Contrast (UPDATE: this is never actually graphed, reserved as comment for now)
    % eval([PFC_type{l},'_',condition,'_Contrast_gamma_struct = ',PFC_type{l},'_',dose_type{1},'_',condition,'_gamma_struct - ',PFC_type{l},'_',dose_type{2},'_',condition,'_gamma_struct;'])

    
    clear PFC gamma_struct Drug_struct Plac_struct
    end

end

%construct Qmoney
for i = 1:length(PFC_type)
    if exist([PFC_type{i},'_Drug_Contrast_gamma_struct'],'var') && ...
       exist([PFC_type{i},'_Plac_Contrast_gamma_struct'],'var')
        eval([PFC_type{i},'_Qmoney_Contrast_gamma_struct = ',PFC_type{i},'_Drug_Contrast_gamma_struct - ',PFC_type{i},'_Plac_Contrast_gamma_struct;'])
    end
end

%construct Drug Green - Plac Green
%note: this will be called Qmoney Green for coding compatability
for i = 1:length(PFC_type)
    if exist([PFC_type{i},'_Drug_Green_gamma_struct'],'var') && ...
       exist([PFC_type{i},'_Plac_Green_gamma_struct'],'var')
        eval([PFC_type{i},'_Qmoney_Green_gamma_struct = ',PFC_type{i},'_Drug_Green_gamma_struct - ',PFC_type{i},'_Plac_Green_gamma_struct;'])
    end
end

%construct Drug Red - Plac Red
%note: this will be called Qmoney Red for coding compatability
for i = 1:length(PFC_type)
    if exist([PFC_type{i},'_Drug_Red_gamma_struct'],'var') && ...
       exist([PFC_type{i},'_Plac_Red_gamma_struct'],'var')
        eval([PFC_type{i},'_Qmoney_Red_gamma_struct = ',PFC_type{i},'_Drug_Red_gamma_struct - ',PFC_type{i},'_Plac_Red_gamma_struct;'])
    end
end

fprintf('done\n')
beep

%% 3) create and save spectrograms

%Note: this function automatically writes in subfolders, beware of this
%      behavior (will not create subfolders if absent)
%Note: uses variables defined in previous cells

% need to manually set frequencies and times if timesout and freqs don't
% exist for ERSP

% need to set subject group, regions, and scale

times = timesout;
frequencies = freqs;

boolUseCustomScaling = true;
c_scale = [-.33 .33];   %dB, lower because 128 channel cap
c_scale_ITC1 = [0 0.3]; %for Green and Red
%c_scale_ITC2 = [0 0.09]; %for Contrasts
c_scale_ITC2 = [-0.04 0.04]; %for Contrasts
%c_scale_ITC3 = [-0.04 0.04]; %for Qmoney only
c_scale_ITC3 = c_scale_ITC2;

%new user defined vars
flagPhase = 0;      %0 = cue,   1 = probe
flagFrequency = 0;  %0 = 4-30Hz,1 = 30-80Hz
strSavePath = '/nfs/erp-modaf/mc/MC-epoched-Feb11-Glenn/Cue/Spectrograms/Test/'; %note: will automatically create subfolders or overwrite
lstConditionsToString = [{'Green'} {'Red'} {'RvG'}];    %used to change naming of condition in filename; comment out to use whatever was in the original lstConditions
%lstConditionsToString = [{'Red'}];
lstConditionsToStringQmoney = [{'DrugGrn-PlacGrn'} {'DrugRed-PlacRed'} {'Qmoney'}]; %only used for Qmoney vars
lstIgnore = [];    %ignore conditions; first dim = dose_type, second dim = condition; each row represents a set, thus cell matrix must be of size n x 2
boolDataIsLog = false;   %enable this if the data was processed with log frequencies (DEBUG_UseLogrithmicFrequencyScaling was set to true)

%gamma_type = [{'Drug_dose3'} {'Plac_single_d3Drug'} {'Contrast'}]; %ms
%gamma_type = [{'Qmoney'}];  %for analysis 2
gamma_type = [{'Plac'} {'Drug'} {'Qmoney'}];    %in any order; fields not found will automatically be skipped
%gamma_type = [{'Plac'}];


%------------------------------------------------------------------------

close all;

%check if data is logrithmic
if boolDataIsLog
    strLogFrequencies = 'native';
else
    strLogFrequencies = 'off';
end
    
for iterator = 1:length(lstConditions)
    condition = lstConditions{iterator};
    
for k = 1:length(gamma_type)
    
    %check if should be ignored
    if ~isempty(lstIgnore)
        idxDoseIgnore = find(cellfun(@(x) strcmp(x,gamma_type{k}),lstIgnore(:,1)));
        idxConditionIgnore = find(cellfun(@(x) strcmp(x,condition),lstIgnore(:,2)));
        boolShouldIgnore = any(intersect(idxDoseIgnore,idxConditionIgnore));
        if boolShouldIgnore
            fprintf('skipped: %s\t%s\n',gamma_type{k},condition)
            continue
        end
    end
    
for j = 1: length(PFC_type)
    
    %check if gamma_type exists as a variable
    if eval(['~exist(''',PFC_type{j},'_',gamma_type{k},'_',condition,'_gamma_struct'',''var'')'])
        fprintf('not found: %s\t%s\n',gamma_type{k},condition)
        break
    end
    
    %Update 1/25/2013: electrode-by-electrode graphing
    if exist('boolUseSingleElectrodes','var') && boolUseSingleElectrodes
        %save original variable somewhere else to keep backwards
        %compatability
        myVar = eval([PFC_type{j},'_',gamma_type{k},'_',condition,'_gamma_struct']);

        %define the position of the electrode, used to find the electrode
        %number
        idxElectrode = 1;

        %extract 1st electrode
        eval([PFC_type{j},'_',gamma_type{k},'_',condition,'_gamma_struct = squeeze(myVar(idxElectrode,:,:));'])
    end
    
    while true %this loop is manually broken at the end of the loop

        %build save path string
        if (flagPhase == 1)
            strSavePath =  strrep(strSavePath,'Cue','Probe');
        end
        strBuilder = strcat(strSavePath,strDataType);
        switch flagFrequency
            case 0
                strBuilder = strcat(strBuilder,' 4-30Hz/');
            case 1
                strBuilder = strcat(strBuilder,' 30-80Hz/');
            otherwise
                error('incorrect flag set');
        end
        switch flagPhase
            case 0
                strBuilder = strcat(strBuilder,'Cue_');
            case 1
                strBuilder = strcat(strBuilder,'Probe_');
            otherwise
                error('incorrect flag set');
        end
        strBuilder = strcat(strBuilder,strDataType,'_');
        switch flagFrequency
            case 0
                strBuilder = strcat(strBuilder,'4-30Hz_');
            case 1
                strBuilder = strcat(strBuilder,'30-80Hz_');
            otherwise
                error('incorrect flag set');
        end
        strBuilder = strcat(strBuilder,strrep(PFC_type{j},'_','-'),'_');
        if strcmp(gamma_type{k},'Qmoney') && exist('lstConditionsToStringQmoney','var')
            strBuilder = strcat(strBuilder,lstConditionsToStringQmoney{iterator});
        elseif ~strcmp(gamma_type{k},'Qmoney') && exist('lstConditionsToString','var')
            strBuilder = strcat(strBuilder,gamma_type{k},'_');
            strBuilder = strcat(strBuilder,lstConditionsToString{iterator});
        else
            strBuilder = strcat(strBuilder,gamma_type{k},'_');
            strBuilder = strcat(strBuilder,lstConditions{iterator});
        end
        if exist('boolUseSingleElectrodes','var') && boolUseSingleElectrodes
            strBuilder = strcat(strBuilder,'_',num2str(eval([PFC_type{j},'(idxElectrode)'])),'.tif');
        else
            strBuilder = strcat(strBuilder,'.tif');
        end

        [strSaveDir strSaveFilename strSaveExtenstion] = fileparts(strBuilder);

        sjset = eval([PFC_type{j},'_',gamma_type{k},'_',condition,'_gamma_struct']);
       
        if strcmp(strDataType,'ITC') && strcmp(gamma_type{k},'Qmoney')
            arrColorScale = c_scale_ITC3;
        elseif strcmp(strDataType,'ITC') && strcmp(condition,'Contrast')
            arrColorScale = c_scale_ITC2;
        elseif strcmp(strDataType,'ITC')
            arrColorScale = c_scale_ITC1;
        else
            arrColorScale = c_scale;
        end
        
        %build title string
        if boolUseCustomScaling
            strTitle = strrep(strSaveFilename,'_','-');
        else
            strTitle = [strrep(strSaveFilename,'_','-'),' (auto scale)'];
        end

        h = tftopo(sjset,timesout,freqs,'title',strTitle,'verbose','off','logfreq',strLogFrequencies);
        if boolUseCustomScaling
            caxis(arrColorScale), colorbar;
        else
            colorbar;
        end
        
        saveas(gcf, strBuilder,'tiff')
        close(gcf);
        
        if ~exist('boolUseSingleElectrodes','var') || ~boolUseSingleElectrodes
            %this is the normal condition with averaged electrodes, exit
            break
        elseif (idxElectrode < size(myVar,1))
            %still more electrodes to process, extract next electrode data
            %and continue
            idxElectrode = idxElectrode +1;
            eval([PFC_type{j},'_',gamma_type{k},'_',condition,'_gamma_struct = squeeze(myVar(idxElectrode,:,:));'])
        else
            %all electrodes have been processed, exit
            break
        end
    end
end
end
end


    