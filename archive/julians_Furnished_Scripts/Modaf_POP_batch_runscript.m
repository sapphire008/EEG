%% Load paths
%% Set Paths
addpath(genpath('/home/gomes/EEGLABv11/'));
addpath(genpath('/nfs/erp-modaf/Furnished Scripts/'))

%% Preprocess cnt files
% Make sure these are "remarked" through Neuroscan
% Make sure broken electrodes are placed
% (raw_eegdir) /nfs/erp-data = Y:\
% set phase to 'delay' or 'probe' 
% set DS to 'DS2' or 'DS1'. usually 'DS2'
% set phase to 'Cue' or 'Probe'
% set dose to 'dose1' ('dose-1' for ms) or 'dose2'or 'dose3'
% place subjects in subj

group = 'ms';
phase = 'Cue'
DS = 'DS2'
dose = 'dose3'
subj = [22 28 30 32]

preprocdir = ['/nfs/erp-modaf/',group,'/',upper(group),'-epoched-Feb11-Glenn/',phase,'/'];
elecs_dir = '/nfs/erp-modaf/elec_files';
raw_eegdir = '/nfs/erp-data/erp2/POP/modaf-re-marked/';
for s = subj
    EEGPreProc_modaf
end


%% run ICA

phase = 'Cue'
dose = 'dose3'
group = 'ms';
preprocdir = ['/nfs/erp-modaf/',group,'/',upper(group),'-epoched-Feb11-Glenn/',phase,'/'];


for s = [22 28 30 32]
runICA_modaf
end

%% View ICs
% set group, dose, and phase
group = 'ms';
dose = 'dose2'
phase = 'Cue'
preprocdir = ['/nfs/erp-modaf/',group,'/',upper(group),'-epoched-Feb11-Glenn/',phase,'/'];

close all

s= [30] 

if s >= 100
        zeroes = ''
    else
        if s >= 10
            zeroes = '0'
        else

            zeroes = '00'
        end
    end

cd(preprocdir);

EEG = pop_loadset( 'filename', [num2str(group),num2str(zeroes),num2str(s),'-',dose,'-PostR-ICA.set'], 'filepath', [preprocdir,num2str(group),num2str(zeroes),num2str(s)]);

% bring up 15 ICs
pop_prop(EEG,0,[1:15]);

% individual components
%eegplot(EEG.icaact,'winlength',[1], 'dispchans', [15]);

% 3-D models
%cd([preprocdir,num2str(group),num2str(zeroes),num2str(s)])
%headplot('setup', EEG.chanlocs, 'headplot.spl')
%pop_headplot(EEG, 0, [1:15], 'Components', [3 5],'load','headplot.spl','view',[180 0],'maplimits','absmax', 'electrodes', 'off')
%% Remove IC's + Artifact reject
% Remember to have created "removed components" and "YG_microsaccades" matrix!

% set group, dose, and preprocdir (from above)
group = 'ms';
phase = 'Cue'
dose = 'dose2'
preprocdir = ['/nfs/erp-modaf/',group,'/',upper(group),'-epoched-Feb11-Glenn/',phase,'/'];

IC_file = [preprocdir,group,'_removed_components_dose.mat']
YG_file = [preprocdir,group,'_YG_microsaccades_dose.mat']


for s = [20]
    ICA_remove_comp_modaf
end


%% Wavelet

addpath('/nfs/erp-modaf/Furnished Scripts/')
clc,clear

c_level = [1 6]; % number of [cycles interval]
freqs = 4:1:30; %frequency bins
times_define = linspace(-200,1200,200);  
%ms: -200,1200,200
%mcpc: -400,1500,200

%set group, phase, and preprocdir (from above)
group = 'ms';
phase = 'Cue';
dose = 'dose1';
preprocdir = ['/nfs/erp-modaf/',group,'/',upper(group),'-epoched-Feb11-Glenn/',phase,'/'];
dirOutput = ['/nfs/erp-modaf/',group,'/',upper(group),'-epoched-Feb11-Glenn/',phase,'/Wavelet results/4-30Hz (no baseline)/'];

%block-wise processing
%Note: this is used for early vs late block analysis. Set the boolean to
%true and list the blocks that should be included in this wavelet
boolUseBlockSeperation = false;
lstBlocks2Process = [3,4];

%subject definitions
sbjs = [2 3 4 5 6 7 8 10 12 13 14 15 16 17 18 20 22 23 24 25 27 28 29 30 32 33 35 36 37 38]; %ms
% sbjs = [2 4 5 6 7 8 9 10 12 14 16 17 18 20 21 22 23 103 108 110 111 115];   %mcpc
% sbjs = [3];   %debug usage
ends = strcat('ms_',phase,'_',dose,'_',strrep(num2str(sbjs(1:length(sbjs)-1:end)),'  ','-'));

%debug parameters (investigational)

DEBUG_UseDS1ForMCPC = false;  %not sure why DS1 is used; data does not seem to suggest that a different set of event codes were used
DEBUG_DisableBaseline = true;   %used for baseline activation investigation in plac red
DEBUG_FLAG_BaselineMethod = 0;  %0 = trial-average, 1 = single-trial, 2 = full; see optional baseline parameters in "help newtimef"
DEBUG_UseLogrithmicFrequencyScaling = false; %used to activate logrithmic scaling of the frequency for scaling of wavelet cycles

%permutation(bootstrapping)
DEBUG_UsePermutation = false;
DEBUG_intAlpha = 0.05;

%--

%autolabel blocks for block-wise processing mode
if boolUseBlockSeperation
    ends = strrep(ends,dose,[dose,'_blk',strrep(num2str(lstBlocks2Process),' ','')]);
end


EEGTimeF_modaf

%% Merge wavelet results
%
%Changelog:
%   3/28/2013:  Added processing of History variable if it exists. If
%               multiple files exist for the same data type (i.e. if
%               subject groups were processed in different batches), then
%               the History from the first file will be used.

clear,clc

%user-defined vars
dirInput = '/nfs/erp-modaf/ms/MS-epoched-Feb11-Glenn/Cue/Wavelet results/4-30Hz (full baseline)/';  %will merge all files with the same prefix (ERSP or ITC) together into new file
dirOutput = '/nfs/erp-modaf/ms/MS-epoched-Feb11-Glenn/Cue/Wavelet results/4-30Hz (full baseline)/';
strFilename = '_ms_dosages_4-30Hz.mat'; %Note: ERSP/ITC will be automatically inserted into beginning of string

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%begin algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Wavelet results merger - MC/PC\nAuthor: Julian Y. Cheng\n\n');

%search input directory
lstInputFiles = dir([dirInput,'*.mat']);
lstInputFiles = struct2cell(lstInputFiles); %convert to struct
lstInputFiles = lstInputFiles(1,:);   %delete unwanted file information

fprintf('Found %i files to process\n\n',length(lstInputFiles))

myERSP = [];
myITC = [];
myFreqs = [];
myTimes = [];
myHistory_ERSP = [];
myHistory_ITC = [];
for cellInputFile = lstInputFiles
    strInputFile = cell2mat(cellInputFile);
    fprintf('Processing %s...',strInputFile)
    
    load([dirInput,strInputFile])
    
    lstStringSplit = regexp(strInputFile,'_','split');
    lstStringSplit(cellfun(@isempty,lstStringSplit)) = [];  %remove empty entries, accounts for double underscores
    
    strDataType = lstStringSplit{1};
    lstVars = who([strDataType,'*']);   %finds the target data variable
    
    %infile integrity check
    if (length(lstVars) ~= 1)   %found more than one instance of data variable, error out
        error('Found more than 1 instance of %s in %s',strDataType,strInputFile)
    end
    if ~exist('freqs','var')
        error('Found no instance of freqs in %s',strInputFile);
    end
    if ~exist('timesout','var')
        error('Found no instance of timesout in %s',strInputFile);
    end
    
    %freqs and timesout integrity check
    if isempty(myFreqs)
        myFreqs = freqs;
    elseif any(myFreqs ~= freqs)
        error('Found incoherent variable freqs in %s',strInputFile);
    end
    if isempty(myTimes)
        myTimes = timesout;
    elseif any(myTimes ~= timesout)
        error('Found incoherent variable timesout in %s',strInputFile);
    end
    
    %first file of each data type, simply store the variable and move on (merge is not performed)
    if isempty(myERSP) && strcmp(strDataType,'ERSP')
        eval(['myERSP = ',lstVars{1},';'])
        if exist('History','var')
            myHistory_ERSP = History;
        end
        fprintf('done\n')
        clearvars -except 'my*' lstInputFiles dirOutput dirInput cellInputFile strFilename
        continue
    elseif isempty(myITC) && strcmp(strDataType,'ITC')
        eval(['myITC = ',lstVars{1},';'])
        if exist('History','var')
            myHistory_ITC = History;
        end
        fprintf('done\n')
        clearvars -except 'my*' lstInputFiles dirOutput dirInput cellInputFile strFilename
        continue
    end
    
    %process file merging
    
    %process filename string
    strPhase = lstStringSplit{3};
    strDose = lstStringSplit{4};
    lstStringSplit(end) = {lstStringSplit{end}(1:end-4)};
    %lstSubjects = cellfun(@(x) str2double(x),lstStringSplit(5:end));    %this parses the filename for subject list, not recommended anymore
    eval(['lstSubjects = find(~cellfun(@isempty,',lstVars{1},'.',strDose,'));'])    %this finds all non-empty entries in the data matrix to find subjects to merge
    
    %merge data
    if strcmp(strDataType,'ERSP')
        eval(['myERSP.',strDose,'(lstSubjects) = ',lstVars{1},'.',strDose,'(lstSubjects);'])
    elseif strcmp(strDataType,'ITC')
        eval(['myITC.',strDose,'(lstSubjects) = ',lstVars{1},'.',strDose,'(lstSubjects);'])
    else
        error('Found unknown data type in %s',strInputFile);
    end
    
    clearvars -except 'my*' lstInputFiles dirOutput dirInput cellInputFile strFilename
    
    fprintf('done\n')
end

fprintf('\nSaving...')

%save variables into compatible name
ERSP = myERSP;
ITC = myITC;
freqs = myFreqs;
timesout = myTimes;

%save mat files
if ~isempty(ERSP)
    fprintf('ERSP...')
    if exist('myHistory_ERSP','var')
        History = myHistory_ERSP;
        save([dirOutput,'ERSP',strFilename],'ERSP','freqs','timesout','History')
    else
        save([dirOutput,'ERSP',strFilename],'ERSP','freqs','timesout')
    end
    clear History
end
if ~isempty(ITC)
    fprintf('ITC...')
    if exist('myHistory_ITC','var')
        History = myHistory_ITC;
        save([dirOutput,'ITC',strFilename],'ITC','freqs','timesout','History')
    else
        save([dirOutput,'ITC',strFilename],'ITC','freqs','timesout')
    end
    clear History
end

fprintf('done\n\nUse this command to check subjects in structures: find(~cellfun(@isempty,[DataType].[Fieldname]))\n\n')

%% Change dosage information to drug/placebo
%
%Changelog:
%   3/28/2013:  Added processing of History variable if it exists. If
%               multiple files exist for the same data type (i.e. if
%               subject groups were processed in different batches), then
%               the History from the first file will be used. Also, changed
%               to UI-based file selection

clear,clc

%user-defined vars
flagStudy = 1;  %0 = MC/PC, 1 = MS dose1/dose2, 2 = MS dose3
dirOutput = '/nfs/erp-modaf/ms/MS-epoched-Feb11-Glenn/Probe/Wavelet results/4-30Hz (full baseline)/';
strOutfilePostfix = '_ms_Drug-Plac_4-30Hz.mat'; %Note: ERSP/ITC will be automatically inserted into beginning of string

%task-specific vars
lstMCPCDose1Drug = [2 4 5 6 7 14 16 17 18 20 23 103 110 112];
lstMSDose1Drug = [4 6 10 12 17 18 20 23 24 29 30 32 33 36];
lstMSDose3Drug = [4 7 10 14 16 17 20 23 24 27 28 29 33 35 36 37 30];
lstMSDose3Plac = [2 8 12 13 15 18 22 25 32 38];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%begin algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%request infile
[strFilename,dirFilepath] = uigetfile('','Select ERSP/ITC wavelet data MAT file');
if isequal(strFilename,0) || isequal(dirFilepath,   0)
    error('file selection canceled, abort')
else
    load([dirFilepath,strFilename])
end

%check if input variable exists
strERSPVar = who('ERSP');
strITCVar = who('ITC');
if xor((length(strERSPVar) ~= 1),(length(strITCVar) ~= 1))
    if (length(strERSPVar) == 1)
        strDataVar = strERSPVar{1};
    else
        strDataVar = strITCVar{1};
    end
else
    error('Input variable error; variable not found or both types co-exist')
end

myData = [];
switch flagStudy
    case 0
        for intSubject = 1:eval(['length(',strDataVar,'.dose1)'])
            if ismember(intSubject,lstMCPCDose1Drug)
                eval(['myData.Drug(intSubject) = ',strDataVar,'.dose1(intSubject);']);
                eval(['myData.Plac(intSubject) = ',strDataVar,'.dose2(intSubject);']);
            else
                eval(['myData.Drug(intSubject) = ',strDataVar,'.dose2(intSubject);']);
                eval(['myData.Plac(intSubject) = ',strDataVar,'.dose1(intSubject);']);
            end
        end
    case 1
        for intSubject = 1:eval(['length(',strDataVar,'.dose1)'])
            if ismember(intSubject,lstMSDose1Drug)
                eval(['myData.Drug(intSubject) = ',strDataVar,'.dose1(intSubject);']);
                eval(['myData.Plac(intSubject) = ',strDataVar,'.dose2(intSubject);']);
            else
                eval(['myData.Drug(intSubject) = ',strDataVar,'.dose2(intSubject);']);
                eval(['myData.Plac(intSubject) = ',strDataVar,'.dose1(intSubject);']);
            end
        end
    case 2
        %note: to preserve subject order in the data vectors, empty cells
        %are inserted
        for intSubject = 1:eval(['length(',strDataVar,'.dose1)'])
            if ismember(intSubject,lstMSDose3Drug)
                eval(['myData.Drug(intSubject) = ',strDataVar,'.dose3(intSubject);']);
                eval('myData.Plac(intSubject) = {[]};');
            else
                eval('myData.Drug(intSubject) = {[]};');
                eval(['myData.Plac(intSubject) = ',strDataVar,'.dose3(intSubject);']);
            end
        end
    otherwise
        error('Invalid flag defined for study')
end

clearvars -except dirOutput strFilename myData strDataVar freqs timesout History strOutfilePostfix

%save file
eval([strDataVar,' = myData;'])
if exist('History','var')
    save([dirOutput,strDataVar,strOutfilePostfix],strDataVar,'freqs','timesout','History')
else
    save([dirOutput,strDataVar,strOutfilePostfix],strDataVar,'freqs','timesout')
end

clear all

beep

%% NaN check
%
%This function is used to check which dimensions have all NaNs; used to
%debug the issue of large regions in the data matrices being NaNs
%Note: this function will list the electrodes that have all NaNs for data
%by subject for each variable

%LOAD IN 1 ERSP/ITC MAT FILE
clc

%user-defined vars
strDose = 'Drug';  %this should be the field in ERSP/ITC; Note: it should not matter which dose is checked, they should all be the same
strCondition = 'incong';    %the condition to test; Note: it should not matter which condition is checked, they should all be the same

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%begin algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get the variable
if exist('ERSP','var')
    myVar = ERSP;
    
    fprintf('Found ERSP\n\n')
elseif exist('ITC','var')
    myVar = ITC;
    
    fprintf('Found ITC\n\n')
else
    error('Cannot find ERSP/ITC variable')
end

%process var of interest
myData = eval(['myVar.',strDose]);
for idxSubject = 1:length(myData)
    if isempty(myData{idxSubject})
        continue
    end
    
    fprintf('Subject %i:\n',idxSubject)

    myElectrodes = eval(['myData{',num2str(idxSubject),'}.',strCondition,'.lead']);
    
    lstNaNs = [];
    lstEmpty = [];
    for idxElectrode = 1:128
        if idxElectrode > length(myElectrodes)
            lstEmpty(end+1) = idxElectrode;
            continue
        elseif isempty(myElectrodes{idxElectrode})
            lstEmpty(end+1) = idxElectrode;
            continue
        end
        
        matrixTemp = isnan(myElectrodes{idxElectrode});
        boolIsAllNaNs = all(matrixTemp(:));

        if boolIsAllNaNs
            lstNaNs(end+1) = idxElectrode;
        end
    end
    
    fprintf('\tEmpty:(%3i) %s\n',length(lstEmpty),strrep(strrep(strrep(num2str(lstEmpty),'    ',','),'   ',','),'  ',','))
    fprintf('\tNaNs: (%3i) %s\n',length(lstNaNs),strrep(strrep(strrep(num2str(lstNaNs),'    ',','),'   ',','),'  ',','))
end

%% Manual baseline-subtraction
%
%This function re-creates baselined data from non-baselined datasets. It is
%only used to investigate why there is activation in the baseline and
%whether the baselining functions in NEWTIMEF() are working correctly.
%
%Note: ITC does not need to be baselined, but code is compatible with
%either data type anyways

%LOAD IN UN-BASELINED ERSP/ITC

%user-defined vars
arrBaseline = [-200 0]; %the baseline range
dirOutput = '/nfs/erp-modaf/mc/MC-epoched-Feb11-Glenn/Cue/Wavelet results/4-30Hz (no baseline subtraction)/';
strFilename = 'mc-pc_Drug-Plac_4-30Hz_baselined.mat';   %automatically appends data type to beginning of string

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%begin algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%parse data type
if exist('ERSP','var') && ~exist('ITC','var')
    strDataVar = 'ERSP';
elseif exist('ITC','var') && ~exist('ERSP','var')
    strDataVar = 'ITC';
else
    error('failed to parse data type; variable not found or both exists')
end

if strcmp(strDataVar,'ITC')
    fprintf('Warning: ITC does not need to be baselined!\n')
end

%find the indexes corresponding to the baseline
idxBaselineStart = find(timesout >= arrBaseline(1),1);
idxBaselineEnd = find(timesout <= arrBaseline(2),1,'last');

fprintf('Processing...\n\n')

%loop through all phases
for cellDose = fieldnames(eval(strDataVar))'
    strDose = cell2mat(cellDose);
        
    %loop through all subjects
    for idxSubject = 1:length(eval([strDataVar,'.',strDose]))
        %skip if subject has no data
        if isempty(eval([strDataVar,'.',strDose,'{idxSubject}']))
            continue
        end
        
        %loop through all conditions
        for cellCondition = fieldnames(eval([strDataVar,'.',strDose,'{idxSubject}']))'
            strCondition = cell2mat(cellCondition);
            
            fprintf('Dose: %s | Subject: %3i | Condition: %-8s | Electrode:',strDose,idxSubject,strCondition)
            
            %loop through electrodes
            %Note: assumes there is a 'lead' fieldname within each
            %condition
            for idxElectrode = 1:length(eval([strDataVar,'.',strDose,'{idxSubject}.',strCondition,'.lead']))
                %skip if electrode has no data
                if isempty(eval([strDataVar,'.',strDose,'{idxSubject}.',strCondition,'.lead{idxElectrode}']))
                    continue
                end
                
                %loop through frequencies
                for idxFrequency = 1:size(eval([strDataVar,'.',strDose,'{idxSubject}.',strCondition,'.lead{idxElectrode}']),1)
                    %get the time series data
                    arrData = eval([strDataVar,'.',strDose,'{idxSubject}.',strCondition,'.lead{idxElectrode}(idxFrequency,:)']);
                    
                    %find the mean of the baseline
                    dblMean = nanmean(arrData(idxBaselineStart:idxBaselineEnd));
                    
                    %subtract the mean from all values
                    arrData = arrData - dblMean;
                    
                    %put the data back into the variable
                    eval([strDataVar,'.',strDose,'{idxSubject}.',strCondition,'.lead{idxElectrode}(idxFrequency,:) = arrData;'])
                end
                
                fprintf('.')
            end
            
            fprintf('\n')
        end
    end
end

fprintf('Saving...')

save([dirOutput,strDataVar,'_',strFilename],strDataVar,'freqs','timesout')

fprintf('done\n')

%% Data extraction by dosage and condition
%
%Note: this function is used to extract the data from the structures (which
%cannot be iterated by indexing) into 4D double arrays. This is to replace
%the current need of doing this individually for every process downstream
%(spectrograms, permutation, headplots, statistics...etc)
%
%Changelog:
%   3/28/2013:  Added processing of History variable if it exists. If
%               multiple files exist for the same data type (i.e. if
%               subject groups were processed in different batches), then
%               the History from the first file will be used.

clear,clc

%user-defined vars
strOutputFolder = '/Processed matrices/';   %will be placed in same folder as input file
strPhase = '';  %leave blank for auto-detect based on input file path, otherwise will use this to form output path
strModifier = '-dose3';   %this is used to modify the "field" label (a.k.a. dosage), used for early vs late analysis and MS dose3
intElectrodes = 128;    %max number of electrodes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%begin algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%request infile
[strFilename,dirFilepath] = uigetfile('','Select ERSP/ITC wavelet data MAT file');
if isequal(strFilename,0) || isequal(dirFilepath,   0)
    error('file selection canceled, abort')
else
    load([dirFilepath,strFilename])
    
    %validate input
    if ~exist('freqs','var') || ~exist('timesout','var')
        error('failed to find variables freqs or timesout')
    end
end

%parse data type
if exist('ERSP','var') && ~exist('ITC','var')
    strDataVar = 'ERSP';
elseif exist('ITC','var') && ~exist('ERSP','var')
    strDataVar = 'ITC';
else
    error('failed to parse data type; variable not found or both exists')
end

%parse phase
if isempty(strPhase)
    if ~isempty(strfind(dirFilepath,'Cue'))
        strPhase = 'Cue';
    elseif ~isempty(strfind(dirFilepath,'Probe'))
        strPhase = 'Phase';
    else
        error('failed to parse phase from path: %s',dirFilepath)
    end
end

%construct output dir and validate
dirOutput = [dirFilepath,strOutputFolder];
if ~exist(dirOutput,'dir')
    fprintf('Warning: output folder does not exist and will be automatically created\n')
    [boolSuccess,strMessage,strMessageID] = mkdir(dirOutput);
    if ~boolSuccess
        error('%s: %s',strMessageID,strMessage)
    end
end
clear boolSuccess strMessage strMessageID

%construct History
if exist('History','var')
    History = fnGetHistory('-except',{'ERSP','ITC','freqs','timesout','History'}, ...
                           'script','Modaf_POP_batch_runscript','source_id',History.id,'source_name',strFilename);
end

%process the data structure
for cellField = eval(['fieldnames(',strDataVar,')'])';  %need to trasponse to interate
    strField = cell2mat(cellField);
    
    %get the data for this field
    dataField = eval([strDataVar,'.',strField]);
    
    %construct list of subjects that have data
    lstSubjects = find(~cellfun(@isempty,dataField));
    
    %remove subjects that are empty
    idxEmptySubjects = find(cellfun(@isempty,dataField));
    dataField(idxEmptySubjects) = [];
    
    %process all subjects
    lstElectrodes = [];
    for idxSubject = 1:length(dataField)
        %create variable template
        myVar = nan(length(dataField),intElectrodes,length(freqs),length(timesout));
        
        %process all sub fields
        for cellSubField = fieldnames(dataField{idxSubject})'
            strSubField = cell2mat(cellSubField);
            
            fprintf('Subject:%4i | Field:%5s | Subfield:%9s | Electrode:',lstSubjects(idxSubject),strField,strSubField)
            
            %get all electrode data
            dataAllElectrodes = eval(['dataField{idxSubject}.',strSubField,'.lead']);
            
            %loop through all electrodes
            dataStorage = nan(size(myVar,2),size(myVar,3),size(myVar,4));
            for idxElectrode = 1:length(dataAllElectrodes)
                %skip if electrode has no data
                if isempty(dataAllElectrodes{idxElectrode})
                    continue
                end
                
                fprintf('.')
                
                %save electrode data
                %if values are complex, convert to magnitude (occurs in ITC)
                if isreal(dataAllElectrodes{idxElectrode})
                    dataStorage(idxElectrode,:,:) = dataAllElectrodes{idxElectrode};
                else
                    dataStorage(idxElectrode,:,:) = abs(dataAllElectrodes{idxElectrode});
                end
            end
            
            fprintf('\n')
            
            %store into subfield-based var
            %note: not the most elegant solution, but since the subfields
            %are contained within the subject dimension there is no
            %alternative
            if ~exist(['data_',strSubField],'var')
                eval(['data_',strSubField,' = myVar;']);
            end
            eval(['data_',strSubField,'(idxSubject,:,:,:) = dataStorage;']);
        end
    end
    
    fprintf('\n')
    
    %find all data storage variables
    lstVars = who('data_*');
    
    %loop through all variables
    for cellVar = lstVars'
        strVar = cell2mat(cellVar);
        
        %parse the subfield name
        strSubField = strrep(strVar,'data_','');
        
        %save variables into structure
        Data.source = strFilename;
        Data.phase = strPhase;
        Data.type = strDataVar;
        Data.field = [strField,strModifier];
        Data.subfield = strSubField;
        Data.subject_list = lstSubjects;
        Data.frequency = freqs;
        Data.time = timesout;
        Data.matrix = eval(strVar);
        
        %save into mat file
        strOutfile = [Data.type,'_',Data.field,'_',Data.subfield,'.mat'];
        if exist([dirOutput,strOutfile],'file')
            fprintf('Saving: %s ... (overwritten) ',strOutfile)
        else
            fprintf('Saving: %s ...',strOutfile)
        end
        save([dirOutput,strOutfile],'Data', 'History')
        fprintf('done\n')
    end
    
    clearvars 'data_*'
end

fprintf('\ndone\n\n')
beep