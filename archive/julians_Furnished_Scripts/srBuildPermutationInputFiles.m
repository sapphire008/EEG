%Generates the necessesary input files for bootstrapping
%Julian Cheng 12/5/2012

%Note: this is a quick fix for the messy MC/PC dataset; not recommended for
%use elsewhere. Copied and modified from SC

%LOAD DESIRED ERSP/ITC BEFORE EXECUTING

%Changelog:
%   1/10/2013:  added functionality to exclude subjects
%   2/5/2013:   removed all abs() calls for ITC; moving this operation to
%               just prior to graphing
%   2/6/2013:   reverted all changes done on 2/5/2013; this is incorrect
%   2/7/2013:   fixed a bug that was causing 0's to be written if electrode
%               is absent in input data matrix (occures only in electrode 128)

%define parameters for extracting data
boolIsERSP = 1;  %1 = ERSP, 0 = ITC in data type
boolIsCue = 1;  %1 = Cue, 0 = Probe for the save path
boolFixFreqsAndTimesout = 0;    %set to 1 if they are missing from input mat (UNTESTED)
strStudy = 'mc-pc';
strFilePostfix = '4-30Hz';
arrFreqRange = [4 30];

% set time_interval and freq_interval in ms
%choose frequency and time range based on ACTUAL ms time
%This will search for it in freqs, timesout
%NOTE: use only if timesout and freqs are available
if ~boolFixFreqsAndTimesout
    freq_range_ms = arrFreqRange;
    time_interval_ms = [timesout(1) timesout(end)];
end

%subject exclusion
boolUseSubjectExclusion = false;
lstExcludedSubjects = [103,108,110,111,115];   %will not process these subjects, use to limit the MC/PC dataset to MC or PC only

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%begin algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Running...');

condition_type = [{'Drug'} {'Plac'}];
trial_type = [{'incong'} {'cong'} {'contrast'}];

%search times and freqs
freq_temp = freqs - freq_range_ms(1);
freq_temp = find(freq_temp>=0);
freq_range(1) = freq_temp(1);
freq_temp = freqs - freq_range_ms(2);
freq_temp = find(freq_temp<=0);
freq_range(2) = freq_temp(end);

times_temp = timesout - time_interval_ms(1);
times_temp = find(times_temp>=0);
time_interval(1) = times_temp(1);
times_temp = timesout - time_interval_ms(2);
times_temp = find(times_temp<=0);
time_interval(2) = times_temp(end);

freq_range = [freq_range(1):freq_range(2)];
time_interval = [time_interval(1):time_interval(2)];

new_freqs = freqs(freq_range);
new_times = timesout(time_interval);

gamma_struct = [];
PFC = 1:128;

if boolIsERSP
    strDataType = 'ERSP';
else
    strDataType = 'ITC';
end

if boolIsCue
    strPhase = 'Cue';
else
    strPhase = 'Probe';
end

if boolFixFreqsAndTimesout
    %get time and freq arrays
    load('/nfs/erp-modaf/mc/MC-epoched-Feb11-Glenn/Cue/old/ITC_mcmc_dose1_4-30Hz_NEW.mat','timesout', 'freqs');
    timesout_185 = timesout;
    freqs_27 = freqs;
    clear timesout freqs
    load('/nfs/erp-modaf/mc/MC-epoched-Feb11-Glenn/Cue/old/ITC_mcmc_dose2_30-80Hz_NEW.mat','timesout', 'freqs');
    timesout_200 = timesout;
    freqs_51 = freqs;
    clear timesout freqs
    
    new_freqs_27 = [4 30];
    new_freqs_51 = [30 80];
    new_freqs_41 = [40 80];
    new_times_185 = timesout_185;
    new_times_200 = timesout_200;

    if eval(['size(',strDataType,'.',condition_type{1},'{1,2}.incong.lead{1},1) == 27'])
        new_freqs = new_freqs_27;
    elseif eval(['size(',strDataType,'.',condition_type{1},'{1,2}.incong.lead{1},1) == 51'])
        new_freqs = new_freqs_51;
    elseif eval(['size(',strDataType,'.',condition_type{1},'{1,2}.incong.lead{1},1) == 41'])
        new_freqs = new_freqs_41;
    else
        error('Unknown freqs')
    end

    if eval(['size(',strDataType,'.',condition_type{1},'{1,2}.incong.lead{1},2) == 185'])
        new_times = new_times_185;
    elseif eval(['size(',strDataType,'.',condition_type{1},'{1,2}.incong.lead{1},2) == 200'])
        new_times = new_times_200;
    else
        error('Unknown timeouts')
    end
end

    
%%get structs

for k= 1:length(condition_type)

    dose = condition_type{k};
    
j=1;
dose_struct_all = eval([strDataType,'.',dose]);
if boolUseSubjectExclusion  %check if subject exclusion is used
    dose_struct_all(lstExcludedSubjects) = [];
    fprintf('Warning: subject exclusion is ON, excluded subjects:\t%s\n',num2str(lstExcludedSubjects))
end
for i=1:length(dose_struct_all)
    if length(dose_struct_all{i}) ~= 0
        eval([dose,'_struct{j} = dose_struct_all{i};'])
        j = j+1;
    end
end

clear i j dose_struct_all


for t = 1:length(trial_type)
    
    title = trial_type{t};

%yield <subject X electrode X freqency X time_interval>

for i=1:length(eval([dose,'_struct']))
     contrast_leads = eval([dose,'_struct{i}.',title,'.lead']);  
    for elec = 1:length(PFC)
        %Update 2/7/2013: check for missing electrodes in data matrix
        %note: sometimes electrode 128 is absent, which will cause
        %erroneous 0's being written
        if elec > length(PFC)
            gamma_struct(i,elec,1:length(freq_range),1:length(time_interval)) = NaN;
            continue
        end
        
        if PFC(elec) <= size(contrast_leads,2) && mean2(contrast_leads{PFC(elec)})~=0 && ~isempty(contrast_leads{PFC(elec)})
            %1. check if this electrode is within the region of interest (in this case though, all regions are of interest becase PFC = 1:128)
            %2. check if the data in this electrode does not average to 0 across frequency and time
            %3. check if this electrode is not empty
            if boolIsERSP
                gamma_struct(i,elec,:,:) = contrast_leads{PFC(elec)}(freq_range,time_interval);
            else 
                gamma_struct(i,elec,:,:) = abs(contrast_leads{PFC(elec)}(freq_range,time_interval));
            end
        else
                gamma_struct(i,elec,1:length(freq_range),1:length(time_interval)) = NaN;
        end
    end
    
end

eval([dose,'_',title,' = gamma_struct;']);

gamma_struct = [];

end

clearvars([dose,'_struct'])

end

%build filename
strFilename = [strDataType,'_',strStudy,'_',strPhase,'_',strFilePostfix,'_TEST.mat'];

clearvars('-except', [condition_type{1},'*'], [condition_type{2},'*'],'name_set','new_freqs','new_times','strFilename','strPhase')

fprintf('Saving...');
save(['/nfs/erp-modaf/mc/MC-epoched-Feb11-Glenn/',strPhase,'/New bootstrapping/',strFilename])
fprintf('done\n');

clear all
beep

%% NaN check
%
%This function is used to check which dimensions have all NaNs; used to
%debug the issue of large regions in the data matrices being NaNs
%Note: this function will list the electrodes that have all NaNs for data
%by subject for each variable

%LOAD IN 1 PERMUTATION INPUT MAT FILE

%user-defined vars
lstSearch = [{'Drug'} {'Plac'}]; %the strings to search for variables

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%begin algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get the variables to search
lstAllVars = who;
lstIndexToDelete = [];
for i = 1:length(lstAllVars)
    boolShouldDelete = true;
    
    for cellSearchString = lstSearch
        strSearch = cell2mat(cellSearchString);
        if ~isempty(strfind(lstAllVars{i},strSearch))
            boolShouldDelete = false;
            break
        end
    end
    
    if boolShouldDelete
        lstIndexToDelete(end+1) = i;
    end
end
lstAllVars(lstIndexToDelete) = [];

%process all vars of interest
for cellVar = lstAllVars'
    strVar = cell2mat(cellVar);
    
    fprintf('Variable: %s\n',strVar)
    
    matrixData = eval(strVar);
    
    for idxSubject = 1:size(matrixData,1)
        fprintf('\tSubject %i:',idxSubject)
        
        intCount = 0;
        for idxElectrode = 1:size(matrixData,2)
            matrixTemp = isnan(squeeze(matrixData(idxSubject,idxElectrode,:,:)));
            boolIsAllNaNs = all(matrixTemp(:));
            
            if boolIsAllNaNs
                fprintf('%i ',idxElectrode)
                intCount = intCount+1;
            end
        end
        
        fprintf('\t(%i)\n',intCount)
    end
    
    fprintf('\n')
end

%% Recreate variables for graphing

%This function is used to debug why the ITC values post permutation are so
%different than the spectrograms. Reformats the data vectors to a format
%that the "Create_Spectrograms_modaf.m" will accept.

%Note: load in permutation input file of interest, then use the "create and
%save spectrograms" code cell in "Create_Spectrograms_modaf.m" to graph.

%user-defined vars
lstSearchTerms = {'Drug*','Plac*'}; %make sure these cover all variables of interest
lstRegions = [{'all_PFC'} {'left_PFC'} {'mid_PFC'} {'right_PFC'} ...
            {'all_Central'} {'left_Central'} {'mid_Central'} {'right_Central'} ...
            {'all_Post'} {'left_Post'} {'mid_Post'} {'right_Post'}];
lstConditions2String = ...  %each pair represents a translation pair, the first string is the target to be replaced, the second string is the replacement
    {'cong','Green'; ...
     'incong','Red'; ...
     'contrast','Contrast'; ...
    };

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%begin algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load in electrode region maps
addpath(genpath('/nfs/erp-modaf/elec_files/'))
load wholehead_elecs.mat
all_PFC = cat(2, wholehead_elecs{1,2},wholehead_elecs{2,2}, wholehead_elecs{3,2});
left_PFC = wholehead_elecs{2,2};
mid_PFC = wholehead_elecs{1,2};
right_PFC = wholehead_elecs{3,2};
all_Central = cat(2, wholehead_elecs{8,2},wholehead_elecs{4,2}, wholehead_elecs{9,2});
left_Central = wholehead_elecs{8,2};
mid_Central = wholehead_elecs{4,2};
right_Central = wholehead_elecs{9,2};
all_Post = cat(2, wholehead_elecs{6,2},wholehead_elecs{5,2}, wholehead_elecs{7,2});
left_Post = wholehead_elecs{6,2};
mid_Post = wholehead_elecs{5,2};
right_Post = wholehead_elecs{7,2};

%search for all variables to process
lstVars = {};
for cellSearch = lstSearchTerms
    strSearch = cell2mat(cellSearch);
    structVarsFound = whos(strSearch);
    if ~isempty(structVarsFound)
        lstVars = [lstVars,{structVarsFound.name}];
    end
end

if isempty(lstVars)
    error('Failed to find any variables to process')
end

%process each variable
for cellVar = lstVars
    strVar = cell2mat(cellVar);
    
    %extract dose and condition information
    strSplit = regexp(strVar,'_','split');
    strDose = strSplit{1};
    strCondition = strSplit{2};
    
    %loop through all regions of interest
    for cellRegion = lstRegions
        strRegion = cell2mat(cellRegion);
        
        %get data in target variable
        dataMaster = eval(strVar);
        
        %get electrode array
        lstElectrodes = eval(strRegion);
        
        %extract subset of data from master
        dataSubset = dataMaster(:,lstElectrodes,:,:);
        
        %average across electrodes then subjects
        dataSubset = squeeze(nanmean(dataSubset,2));
        dataSubset = squeeze(nanmean(dataSubset,1));
        
        %build variable name
        strBuilder = [strRegion,'_',strDose,'_'];
        idxCondition2String = find(cellfun(@(x) strcmp(x,strCondition),lstConditions2String(:,1)));   %find the condition string in the first dimension of lstConditions2String
        if isempty(idxCondition2String)
            fprintf('Warning: failed to find entry in condition to string list, using original name for: %s %s %s\n',strRegion,strDose,strCondition);
            strBuilder = [strBuilder,strCondition,'_gamma_struct'];
        else
            strBuilder = [strBuilder,lstConditions2String{idxCondition2String,2},'_gamma_struct'];
        end
        
        %create the variable
        eval([strBuilder,' = dataSubset;'])
    end
end

%construct Qmoney
for i = 1:length(lstRegions)
    if exist([lstRegions{i},'_Drug_Contrast_gamma_struct'],'var') && ...
       exist([lstRegions{i},'_Plac_Contrast_gamma_struct'],'var')
        eval([lstRegions{i},'_Qmoney_Contrast_gamma_struct = ',lstRegions{i},'_Drug_Contrast_gamma_struct - ',lstRegions{i},'_Plac_Contrast_gamma_struct;'])
    end
end

%construct Drug Green - Plac Green
%note: this will be called Qmoney Green for coding compatability
for i = 1:length(lstRegions)
    if exist([lstRegions{i},'_Drug_Green_gamma_struct'],'var') && ...
       exist([lstRegions{i},'_Plac_Green_gamma_struct'],'var')
        eval([lstRegions{i},'_Qmoney_Green_gamma_struct = ',lstRegions{i},'_Drug_Green_gamma_struct - ',lstRegions{i},'_Plac_Green_gamma_struct;'])
    end
end

%construct Drug Red - Plac Red
%note: this will be called Qmoney Red for coding compatability
for i = 1:length(lstRegions)
    if exist([lstRegions{i},'_Drug_Red_gamma_struct'],'var') && ...
       exist([lstRegions{i},'_Plac_Red_gamma_struct'],'var')
        eval([lstRegions{i},'_Qmoney_Red_gamma_struct = ',lstRegions{i},'_Drug_Red_gamma_struct - ',lstRegions{i},'_Plac_Red_gamma_struct;'])
    end
end

%set variables for graphing code compatability

if ~isempty(strfind(strFilename,'ERSP'))
    strDataType = 'ERSP';
elseif ~isempty(strfind(strFilename,'ITC'))
    strDataType = 'ITC';
else
    error('Could not define datatype')
end

timesout = new_times;
freqs = new_freqs;
lstConditions = lstConditions2String(:,2);
PFC_type = lstRegions;

beep