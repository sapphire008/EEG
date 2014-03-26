%% Data processing for headplot generation
%by: Julian Y. Cheng    12/19/2012
%
%Adapted from SC script "Headmap.m" for MC/PC
%
%This script conditions the ERSP outputs from the wavelet analysis to
%produce datasets that can be plotted directly or can undergo permutation
%(bootstrapping) before plotting.

%% Extract electrode data from ERSP
%Produces <subject X electrode> matrix
%
%MUST LOAD IN 1 COPY OF DESIRED ERSP

clc

%user-defined vars
strPhase = 'Cue';
lstFrequencyRange = [15 30];  %Hz
lstTimeInterval = [0 250];  %ms

%subject exclusion
boolUseSubjectExclusion = true;
lstExcludedSubjects = [103,108,110,111,115];   %will not process these subjects, use to limit the MC/PC dataset to MC or PC only

%task-specific vars
intElectrodes = 128;
lstConditions = [{'Drug'} {'Plac'}];
lstTrialTypes = [{'incong'} {'cong'} {'contrast'}];
dirOutput = ['/nfs/erp-modaf/mc/MC-epoched-Feb11-Glenn/',strPhase,'/Headplots/'];
strFilename = ['HDPLT_ERSP_MC-only_',strPhase,'_',num2str(lstFrequencyRange(1)),'-',num2str(lstFrequencyRange(2)),'Hz_',num2str(lstTimeInterval(1)),'-',num2str(lstTimeInterval(2)),'ms.mat'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%begin algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Data processing for headplot generation - MC/PC\nAuthor: Julian Y. Cheng\n\n');

%search times and freqs
freq_temp = freqs - lstFrequencyRange(1);
freq_temp = find(freq_temp>=0);
freq_range(1) = freq_temp(1);
freq_temp = freqs - lstFrequencyRange(2);
freq_temp = find(freq_temp<=0);
freq_range(2) = freq_temp(end);
times_temp = timesout - lstTimeInterval(1);
times_temp = find(times_temp>=0);
time_interval(1) = times_temp(1);
times_temp = timesout - lstTimeInterval(2);
times_temp = find(times_temp<=0);
time_interval(2) = times_temp(end);
arrFrequencyRange = freq_range(1):freq_range(2);
arrTimeInterval = time_interval(1):time_interval(2);

%process ERSP
fprintf('Processing...\n')
for i = 1:length(lstConditions)
    strCondition = lstConditions{i};
    fprintf('\t%s...\n',strCondition);
    
    for j = 1:length(lstTrialTypes)
        strTrialType = lstTrialTypes{j};
        fprintf('\t\t%s...',strTrialType)
        
        %get all subject electrode data from ERSP
        eval(['structAllSubjects = ERSP.',strCondition,';'])
        if boolUseSubjectExclusion  %check if subject exclusion is used
            structAllSubjects(lstExcludedSubjects) = [];
            fprintf('Warning: subject exclusion is ON, excluded subjects:\t%s\n',num2str(lstExcludedSubjects))
        end
        idxEmptyEntries = cellfun(@isempty,structAllSubjects);
        structAllSubjects(idxEmptyEntries) = [];
        
        %process into <subject X electrode X time interval> matrix
        matrixData = NaN(length(structAllSubjects),intElectrodes,length(arrTimeInterval));
        for idxSubject = 1:size(matrixData,1)
            eval(['dataElectrodes = structAllSubjects{',num2str(idxSubject),'}.',strTrialType,'.lead;'])
            
            if (length(dataElectrodes) ~= intElectrodes)
                %report subjects that have fewer than maximum electrodes
                %Note: this should not happen, because even excluded
                %electrodes are included in data matrix as empty entries to
                %preserve electrode position. Not sure why this happens.
                fprintf('%i(%i)..',idxSubject,length(dataElectrodes))
            end
            
            for idxElectrode = 1:length(dataElectrodes)
                eval(['dataThisElectrode = dataElectrodes{',num2str(idxElectrode),'};'])  %frequency x time data for a single electrode
                
                if isempty(dataThisElectrode),continue,end
                
                %average accross frequency
                matrixData(idxSubject,idxElectrode,1:length(arrTimeInterval)) = nanmean(dataThisElectrode(arrFrequencyRange,arrTimeInterval),1);
            end
        end
        
        %average across time
        matrixData = squeeze(nanmean(matrixData,3));    %yields <subject X electrode>
        
        eval(['matrixSubjectByElectrode_',strCondition,'_',strTrialType,' = matrixData;']);
        
        fprintf('done\n')
    end
end

%construct Qmoney
if exist('matrixSubjectByElectrode_Drug_contrast','var') && ...
   exist('matrixSubjectByElectrode_Plac_contrast','var')
    matrixSubjectByElectrode_Qmoney_contrast = matrixSubjectByElectrode_Drug_contrast - matrixSubjectByElectrode_Plac_contrast;
end

%Drug Green - Plac Green not implemented

fprintf('done\n\nSaving...')

clearvars -except 'matrixSubjectByElectrode*' dirOutput strFilename strPhase

save([dirOutput,strFilename])

fprintf('done\n')

clear all

%% Plot topographic map

%LOAD DATA MATRIX/SIG MATRIX FIRST

%user-defined vars
boolIsSigMatrix = 0;             %0 = 'data matrix', 1 = 'sig matrix'; changes how the data is adjusted before graphing
intRejectLimit = 17;            %maximum number of subjects that can be excluded per electrode; if exceeded the electrode will be excluded. Not used if boolIsSigMatrix is set
%strVarSearch = '*Plac_contrast';    %a common part of the input variable used to find the proper input data; not used if boolIsSigMatrix is set
strVarSearch = '*Qmoney_contrast';  %MAKE SURE THIS MATCHES WITH VARIABLE strFilePrefix
%strVarSearch = '*Plac_cong';

%autoname
boolUseAutoname = 1;    %automatically saves figures and closes them, format: [strFilePrefix]_[flagFrequency]_[strTimeInterval].tif
dirOutput = '/nfs/erp-modaf/mc/MC-epoched-Feb11-Glenn/Cue/Headplots/Figures/';
%strFilePrefix = 'ERSP_Cue_Plac_RvG';
strFilePrefix = 'ERSP_Cue_Qmoney';
%strFilePrefix = 'ERSP_Cue_Plac_Green';
flagFrequency = 2;      %0 = 4-8Hz, 1 = 8-15Hz, 2 = 15-30Hz
strTimeInterval = '0-250ms';

%task-specific vars
pathCleanEEGSet = '/nfs/erp-modaf/mc/ERP/Generic_set_1trial.set';   %needs to have complete chanloc information for all channels

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%begin algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get chanlocs from a complete .set file
structEEG = pop_loadset(pathCleanEEGSet);
structChanlocs = structEEG.chanlocs;
    
%condition data (remove bad electrodes)
idxElectrodes = 1:128;
if boolIsSigMatrix
    %find electrodes to exclude
    idxElectrodes2Remove = find(isnan(actual_t_matrix));
    
    dataElectrodes = sig_matrix;
else
    %find the input variable
    lstInputVars = who(strVarSearch);
    
    %we only expect one, error out if found more than 1
    if (length(lstInputVars) > 1)
        error('more than 1 input variable found; see lstInputVars')
    elseif isempty(lstInputVars)
        error('no variables found');
    else
        strInputVar = cell2mat(lstInputVars);
        eval(['dataElectrodes = ',strInputVar,';']);
    end
    
    %find electrodes to exclude
    idxElectrodes2Remove = [];
    for i = 1:length(dataElectrodes)
        intExcludedSubjects = sum(isnan(dataElectrodes(:,i)));
        if (intExcludedSubjects > intRejectLimit)
            idxElectrodes2Remove(end+1) = i;
            fprintf('Removed electrode %i for having %i excluded subjects\n',i,intExcludedSubjects);
        end
    end
    
    dataElectrodes = nanmean(dataElectrodes,1);
end
idxElectrodes(idxElectrodes2Remove) = [];

figure,colorbar;
hdlePlot = topoplot(dataElectrodes,structChanlocs,'electrodes','on','emarker',{'.','k',5,1},'plotchans',idxElectrodes);

if boolUseAutoname
    %build output file path
    switch flagFrequency
        case 0
            strFrequency = '4-8Hz';
        case 1
            strFrequency = '8-15Hz';
        case 2
            strFrequency = '15-30Hz';
        otherwise
            error('Invalid value for flagFrequency')
    end
    strFilename = [strFilePrefix,'_',strFrequency,'_',strTimeInterval,'.tif'];
    
    %save figure and close
    saveas(gcf,[dirOutput,strFilename],'tiff')
    close(gcf)
    beep
end

clear all

%% Multiplot for sig matrix

clear all
close all
clc

%user-defined vars
dirInput = '/nfs/erp-modaf/mc/MC-epoched-Feb11-Glenn/Cue/Headplots/Bootstrap/';   %will process all files contained within, excluding subfolders
dirOutput = '/nfs/erp-modaf/mc/MC-epoched-Feb11-Glenn/Cue/Headplots/Figures/';

%task-specific vars
pathCleanEEGSet = '/nfs/erp-modaf/mc/ERP/Generic_set_1trial.set';   %needs to have complete chanloc information for all channels

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%begin algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

warning off backtrace

%search input directory
lstInputFiles = dir([dirInput,'*.mat']);
lstInputFiles = struct2cell(lstInputFiles); %convert to struct
lstInputFiles = lstInputFiles(1,:);   %delete unwanted file information

%get chanlocs from a complete .set file
structEEG = pop_loadset(pathCleanEEGSet);
structChanlocs = structEEG.chanlocs;

for cellInputFile = lstInputFiles
    strInputFile = cell2mat(cellInputFile);
    load([dirInput,strInputFile])
    
    %check infile integrity
    if ~exist('sig_matrix','var')
        warning('Cannot find sig_matrix: %s',strInputFile)
        continue
    end
    
    %condition data (remove bad electrodes)
    
    %find electrodes to exclude
    idxElectrodes = 1:128;
    idxElectrodes2Remove = find(isnan(actual_t_matrix));

    dataElectrodes = sig_matrix;
    idxElectrodes(idxElectrodes2Remove) = [];
    
    %plot figure and save

    figure,colorbar;
    hdlePlot = topoplot(dataElectrodes,structChanlocs,'electrodes','on','emarker',{'.','k',5,1},'plotchans',idxElectrodes);

    %build save path
    
    %parse source file string
    lstStringSplit = regexp(strSourceFile,'_','split');
    strDataType = lstStringSplit{2};
    strPhase = lstStringSplit{4};
    strFrequency = lstStringSplit{5};
    strTimeInterval = lstStringSplit{6}(1:end-4);
    
    %parse infile string
    lstStringSplit = regexp(strInputFile,'_','split');
    if strcmp(lstStringSplit{5},'Qmoney')    %exception for Qmoney
        strCondition = [lstStringSplit{5}];
    else
        strCondition = [lstStringSplit{5},'_',lstStringSplit{6}];
    end
    strFilePostfix = cell2mat(regexp(lstStringSplit{end},'\d+','match'));
    if ~isempty(strFilePostfix)
        strFilePostfix = ['_',strFilePostfix];
    end
    
    pathOutfile = [dirOutput,strDataType,'_',strPhase,'_',strCondition,'_',strFrequency,'_',strTimeInterval,'_bootstrap',strFilePostfix,'.tif'];
        
    saveas(gcf,pathOutfile,'tiff')
    close(gcf)
end

clear all

fprintf('\ndone\n')
beep

warning on backtrace

%% Plot 3-D head map

%LOAD DATA MATRIX/SIG MATRIX FIRST

%user-defined vars
boolIsSigMatrix = 0;             %0 = 'data matrix', 1 = 'sig matrix'; changes how the data is adjusted before graphing
intRejectLimit = 17;            %maximum number of subjects that can be excluded per electrode; if exceeded the electrode will be excluded. Not used if boolIsSigMatrix is set
%strVarSearch = '*Plac_contrast';    %a common part of the input variable used to find the proper input data; not used if boolIsSigMatrix is set
strVarSearch = '*Qmoney_contrast';  %MAKE SURE THIS MATCHES WITH VARIABLE strFilePrefix
%strVarSearch = '*Plac_cong';

%autoname
boolUseAutoname = 0;    %automatically saves figures and closes them, format: [strFilePrefix]_[flagFrequency]_[strTimeInterval].tif
dirOutput = '/nfs/erp-modaf/mc/MC-epoched-Feb11-Glenn/Cue/Headplots/Figures/';
%strFilePrefix = 'ERSP_Cue_Plac_RvG';
strFilePrefix = 'ERSP_Cue_Qmoney';
%strFilePrefix = 'ERSP_Cue_Plac_Green';
flagFrequency = 2;      %0 = 4-8Hz, 1 = 8-15Hz, 2 = 15-30Hz
strTimeInterval = '0-250ms';

%task-specific vars
pathCleanEEGSet = '/nfs/erp-modaf/mc/ERP/Generic_set_1trial.set';   %needs to have complete chanloc information for all channels
pathSpline = '/nfs/erp-modaf/mc/MC-epoched-Feb11-Glenn/Cue/Headplots/HDPLT_SPLINE.spl';  %spline file will be saved here

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%begin algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get chanlocs from a complete .set file
structEEG = pop_loadset(pathCleanEEGSet);
    
%condition data (remove bad electrodes)
if boolIsSigMatrix
    %find electrodes to exclude
    idxElectrodes2Remove = find(isnan(actual_t_matrix));
    
    dataElectrodes = sig_matrix;
else
    %find the input variable
    lstInputVars = who(strVarSearch);
    
    %we only expect one, error out if found more than 1
    if (length(lstInputVars) > 1)
        error('more than 1 input variable found; see lstInputVars')
    elseif isempty(lstInputVars)
        error('no variables found');
    else
        strInputVar = cell2mat(lstInputVars);
        eval(['dataElectrodes = ',strInputVar,';']);
    end
    
    %find electrodes to exclude
    idxElectrodes2Remove = [];
    for i = 1:length(dataElectrodes)
        intExcludedSubjects = sum(isnan(dataElectrodes(:,i)));
        if (intExcludedSubjects > intRejectLimit)
            idxElectrodes2Remove(end+1) = i;
            fprintf('Removed electrode %i for having %i excluded subjects\n',i,intExcludedSubjects);
        end
    end
    
    dataElectrodes = nanmean(dataElectrodes,1);
end
structEEG = pop_select(structEEG,'nochannel',idxElectrodes2Remove);
dataElectrodes(idxElectrodes2Remove) = [];

%create spline file
headplot('setup',structEEG.chanlocs,pathSpline)

hdlePlot = headplot(dataElectrodes,pathSpline,'view',[0 90],'cbar',0,'electrodes','on');

if boolUseAutoname
    %build output file path
    switch flagFrequency
        case 0
            strFrequency = '4-8Hz';
        case 1
            strFrequency = '8-15Hz';
        case 2
            strFrequency = '15-30Hz';
        otherwise
            error('Invalid value for flagFrequency')
    end
    strFilename = [strFilePrefix,'_',strFrequency,'_',strTimeInterval,'.tif'];
    
%     %save figure and close
%     saveas(gcf,[dirOutput,strFilename],'tiff')
%     close(gcf)
%     beep
end

%% Plot (MUST re-run previous step EACH time!) OLD, for 3-D headplots

% set scale
scale_limits = [-1.5 1.5];

EEG = pop_loadset('/nfs/erp-stroop/Glenn-Mar11/Cued_Stroop/sc-generic_set.set');

% remove electrode locations by matching electrode number with EEG.chanlocs
%.type

bad_elecs = find(isnan(cue_mean_gamma_struct(:,1))');

elecs = cellfun(@str2num,{EEG.chanlocs.type});

remove_elecs = [];
for b=1:length(bad_elecs)
    if ~isempty(find(elecs == bad_elecs(b)));
        remove_elecs(end+1) = find(elecs == bad_elecs(b));
    end
end

EEG = pop_select( EEG, 'nochannel', remove_elecs);

cue_mean_gamma_struct(bad_elecs,:) = [];

EEG.data = mean(cue_mean_gamma_struct,2);

headplot('setup', EEG.chanlocs, 'headplot.spl')

figure; headplot(EEG.data,'headplot.spl','view',[0 90],'maplimits',scale_limits,'cbar',0,'electrodes','off')    

