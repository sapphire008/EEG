%Permutation test pipeline
%Author: Julian Y. Cheng
%2/6/2013
%
%This script is an adaption of the srBuildPermutationInputFiles and
%srPermutation in an attempt to clean up the process. This script uses
%processed inputs, not the raw wavelet output
%
%Note: this script makes use of the new permutation_test function by Scott;
%it requires that this script be run on a machine equiped with a NVIDIA GPU
%for CUDA
%
%Changelog:
%   3/28/2013:  Added processing of History variable if it exists. If
%               multiple files exist for the same data type (i.e. if
%               subject groups were processed in different batches), then
%               the History from the first file will be used.

clc,clear

%add path of GPU function
addpath /nfs/pkg64/gpu/

%user-defined vars---------------------------------------------------------

%general
flagDataType = 1;   %1 = ERSP, 2 = ITC
flagPhase = 1;      %1 = Cue, 2 = probe
flagTTest = 2;      %1 = 1-tailed, 2 = 2-tailed
strInputFolder = '/Wavelet results/4-30Hz (full baseline)/Processed matrices/'; %combined with dirRoot and flagPhase to make the input search path
strOutputFolder = '/Permutation/ERSP 4-30Hz paired 2-tail (full baseline)/';    %combined with dirRoot and flagPhase to make save path

%permutation parameters
intCycles = 4000000;  %the number of permutation cycles
dblAlpha = 0.05;    %the alpha value to use in selecting p-values

%comparisons to run
%note: each entry in the structure array is a comparison group with the
%following fieldnames:
%   .name = the name of the comparison, will be used in save filename
%   .pos.dose = the dosage of the first group
%   .pos.cond = the condition of the first group
%   .neg.dose = the dosage of the second group
%   .neg.cond = the condition of the secound group
%*dosages and conditions must match labels used in processed wavelet files
%*comparisons with required files missing will be automatically ignored
%[all blocks]
if true
    structComparisons(1).name = 'Drug_RvG';
    structComparisons(1).pos.dose = 'Drug';
    structComparisons(1).pos.cond = 'incong';
    structComparisons(1).neg.dose = 'Drug';
    structComparisons(1).neg.cond = 'cong';
    structComparisons(2).name = 'Plac_RvG';
    structComparisons(2).pos.dose = 'Plac';
    structComparisons(2).pos.cond = 'incong';
    structComparisons(2).neg.dose = 'Plac';
    structComparisons(2).neg.cond = 'cong';
    structComparisons(3).name = 'Qmoney';
    structComparisons(3).pos.dose = 'Drug';
    structComparisons(3).pos.cond = 'contrast';
    structComparisons(3).neg.dose = 'Plac';
    structComparisons(3).neg.cond = 'contrast';
end
%[early blocks]
if true
    structComparisons(4).name = 'Drug_RvG_early';
    structComparisons(4).pos.dose = 'Drug-early';
    structComparisons(4).pos.cond = 'incong';
    structComparisons(4).neg.dose = 'Drug-early';
    structComparisons(4).neg.cond = 'cong';
    structComparisons(5).name = 'Plac_RvG_early';
    structComparisons(5).pos.dose = 'Plac-early';
    structComparisons(5).pos.cond = 'incong';
    structComparisons(5).neg.dose = 'Plac-early';
    structComparisons(5).neg.cond = 'cong';
    structComparisons(6).name = 'Qmoney_early';
    structComparisons(6).pos.dose = 'Drug-early';
    structComparisons(6).pos.cond = 'contrast';
    structComparisons(6).neg.dose = 'Plac-early';
    structComparisons(6).neg.cond = 'contrast';
end
%[late blocks]
if true
    structComparisons(7).name = 'Drug_RvG_late';
    structComparisons(7).pos.dose = 'Drug-late';
    structComparisons(7).pos.cond = 'incong';
    structComparisons(7).neg.dose = 'Drug-late';
    structComparisons(7).neg.cond = 'cong';
    structComparisons(8).name = 'Plac_RvG_late';
    structComparisons(8).pos.dose = 'Plac-late';
    structComparisons(8).pos.cond = 'incong';
    structComparisons(8).neg.dose = 'Plac-late';
    structComparisons(8).neg.cond = 'cong';
    structComparisons(9).name = 'Qmoney_late';
    structComparisons(9).pos.dose = 'Drug-late';
    structComparisons(9).pos.cond = 'contrast';
    structComparisons(9).neg.dose = 'Plac-late';
    structComparisons(9).neg.cond = 'contrast';
end
%[blocks 1,2]
if true
    structComparisons(10).name = 'Drug_RvG_blk12';
    structComparisons(10).pos.dose = 'Drug-blk12';
    structComparisons(10).pos.cond = 'incong';
    structComparisons(10).neg.dose = 'Drug-blk12';
    structComparisons(10).neg.cond = 'cong';
    structComparisons(11).name = 'Plac_RvG_blk12';
    structComparisons(11).pos.dose = 'Plac-blk12';
    structComparisons(11).pos.cond = 'incong';
    structComparisons(11).neg.dose = 'Plac-blk12';
    structComparisons(11).neg.cond = 'cong';
    structComparisons(12).name = 'Qmoney_blk12';
    structComparisons(12).pos.dose = 'Drug-blk12';
    structComparisons(12).pos.cond = 'contrast';
    structComparisons(12).neg.dose = 'Plac-blk12';
    structComparisons(12).neg.cond = 'contrast';
end
%[blocks 3,4]
if true
    structComparisons(13).name = 'Drug_RvG_blk34';
    structComparisons(13).pos.dose = 'Drug-blk34';
    structComparisons(13).pos.cond = 'incong';
    structComparisons(13).neg.dose = 'Drug-blk34';
    structComparisons(13).neg.cond = 'cong';
    structComparisons(14).name = 'Plac_RvG_blk34';
    structComparisons(14).pos.dose = 'Plac-blk34';
    structComparisons(14).pos.cond = 'incong';
    structComparisons(14).neg.dose = 'Plac-blk34';
    structComparisons(14).neg.cond = 'cong';
    structComparisons(15).name = 'Qmoney_blk34';
    structComparisons(15).pos.dose = 'Drug-blk34';
    structComparisons(15).pos.cond = 'contrast';
    structComparisons(15).neg.dose = 'Plac-blk34';
    structComparisons(15).neg.cond = 'contrast';
end
%[blocks 5,6]
if true
    structComparisons(16).name = 'Drug_RvG_blk56';
    structComparisons(16).pos.dose = 'Drug-blk56';
    structComparisons(16).pos.cond = 'incong';
    structComparisons(16).neg.dose = 'Drug-blk56';
    structComparisons(16).neg.cond = 'cong';
    structComparisons(17).name = 'Plac_RvG_blk56';
    structComparisons(17).pos.dose = 'Plac-blk56';
    structComparisons(17).pos.cond = 'incong';
    structComparisons(17).neg.dose = 'Plac-blk56';
    structComparisons(17).neg.cond = 'cong';
    structComparisons(18).name = 'Qmoney_blk56';
    structComparisons(18).pos.dose = 'Drug-blk56';
    structComparisons(18).pos.cond = 'contrast';
    structComparisons(18).neg.dose = 'Plac-blk56';
    structComparisons(18).neg.cond = 'contrast';
end
%[blocks 7,8]
if true
    structComparisons(19).name = 'Drug_RvG_blk78';
    structComparisons(19).pos.dose = 'Drug-blk78';
    structComparisons(19).pos.cond = 'incong';
    structComparisons(19).neg.dose = 'Drug-blk78';
    structComparisons(19).neg.cond = 'cong';
    structComparisons(20).name = 'Plac_RvG_blk78';
    structComparisons(20).pos.dose = 'Plac-blk78';
    structComparisons(20).pos.cond = 'incong';
    structComparisons(20).neg.dose = 'Plac-blk78';
    structComparisons(20).neg.cond = 'cong';
    structComparisons(21).name = 'Qmoney_blk78';
    structComparisons(21).pos.dose = 'Drug-blk78';
    structComparisons(21).pos.cond = 'contrast';
    structComparisons(21).neg.dose = 'Plac-blk78';
    structComparisons(21).neg.cond = 'contrast';
end

%regions to run
%note: each region has an accompanying switch, turn them on/off to select
%the desired regions to be run
switchRegions = { ...
    'all-PFC'       '0'; ...
    'left-PFC'      '1'; ...
    'mid-PFC'       '1'; ...
    'right-PFC'     '1'; ...
    'all-Central'   '0'; ...
    'left-Central'  '1'; ...
    'mid-Central'   '1'; ...
    'right-Central' '1'; ...
    'all-Post'      '0'; ...
    'left-Post'     '0'; ...
    'mid-Post'      '0'; ...
    'right-Post'    '0'; ...
    };

%frequency settings
boolAutoDetect_Frequency = true;    %auto-detect is based on the min and max values of the Data.frequency variable; if disabled use arrFrequency to define range
arrFrequency = [4 30];              %ignored if auto-detect is enabled

%time scale settings
boolAutoDetect_Time = true;    %auto-detect is based on the min and max values of the Data.times variable; if disabled use arrTimes to define range
arrTimes = [-400 1500];

%subject exclusion
%note: this is used to limit processing to a subset of the subject pool,
%for instance mc-only in the mc-pc dataset
boolUseSubjectExclusion = false;
lstExcludedSubjects = [103,107,110,111,115];

%task-specific vars--------------------------------------------------------

%general
dirRoot = '/nfs/erp-modaf/ms/MS-epoched-Feb11-Glenn/';
strStudy = 'ms';                                                      %used in save file name
strMyName = 'pipePermutationTest';                                      %name of this file, used in history
strFilePrefix = 'PRMT';                                                 %the prefix added to the beginning of the output filename
pathElectrodeMap = '/nfs/erp-modaf/elec_files/wholehead_elecs.mat';     %defines the regions
intElectrodes = 128;                                                    %max number of electrodes in data

%legacy function call
%note: param1 = row vector of positive term
%      param2 = row vector of negative term
%      param3 = test type
%      param4 = permutation cycles
%      output1 = p-value (range[0,1])
%      output2 = t-value
strFunction = 'permutation_test';   %function name
strTestType = 'paired';             %parameter used to choose permutation method; paired = paired t-test (non-bootstrapping), pooled = mean-difference bootstrapping

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%begin algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Permutation test pipeline: %s\nAuthor: Julian Y. Cheng\n\n',strStudy);

%parse data type
switch flagDataType
    case 1  %ERSP
        strDataType = 'ERSP';
    case 2  %ITC
        strDataType = 'ITC';
    otherwise
        error('incorrect flag set for data type')
end
clear flagDataType

%parse phase
switch flagPhase
    case 1  %cue
        strPhase = 'Cue';
    case 2  %probe
        strPhase = 'Probe';
    otherwise
        error('incorrect flag set for phase')
end
clear flagPhase

%parse t test tail
switch flagTTest
    case 1  %1-tailed
        boolUseSingleTail = true;
    case 2  %2-tailed
        boolUseSingleTail = false;
    otherwise
        error('incorrect flag set for t test tail')
end
clear flagTTest

%find processed wavelet data (only load 1st one for parsing parameters, all data loaded after history is saved)
dirInfolder = [dirRoot,strPhase,strInputFolder];
if ~exist(dirInfolder,'dir')
    error('failed to find input folder: %s',dirInfolder)
end
lstInfiles = dir([dirInfolder,strDataType,'*']);
lstInfiles = {lstInfiles.name};
load([dirInfolder,lstInfiles{1}])

%parse frequency and find indexes(keep)
if boolAutoDetect_Frequency
    arrFrequency = [min(Data.frequency),max(Data.frequency)];
end
arrFrequencyIndexes = find(Data.frequency >= arrFrequency(1),1):find(Data.frequency <= arrFrequency(2),1,'last');

%parse timepoints and find indexes(keep)
if boolAutoDetect_Time
    arrTimes = [min(Data.time),max(Data.time)];
end
arrTimesIndexes = find(Data.time >= arrTimes(1),1):find(Data.time <= arrTimes(2),1,'last');

clear Data History

%generate output directory (validation after confirm user input)
dirOutput = [dirRoot,strPhase,strOutputFolder];

%cleanup comparisons to only the ones that have required files found
lstDeleteIndexes = [];
for i = 1:length(structComparisons)
    %exclude this comparison if any required file is not found
    if all(cellfun(@isempty,regexp(lstInfiles,[structComparisons(i).pos.dose,'_',structComparisons(i).pos.cond])))
        lstDeleteIndexes(end+1) = i;
    elseif all(cellfun(@isempty,regexp(lstInfiles,[structComparisons(i).neg.dose,'_',structComparisons(i).neg.cond])))
        lstDeleteIndexes(end+1) = i;
    end
end
structComparisons(lstDeleteIndexes) = [];
clear lstDeleteIndexes i

%construct history
myHistory = fnGetHistory('script',strMyName);

%info messages
fprintf('\n[Parameters]:\n')
fprintf('Input dir: %s\n',dirInfolder)
fprintf('Output dir: %s\n',dirOutput)
fprintf('Phase: %s\n',strPhase)
fprintf('Data type: %s\n',strDataType)
fprintf('Frequency: %i-%iHz\n',arrFrequency(1),arrFrequency(2))
fprintf('Time interval: %i-%ims\n',arrTimes(1),arrTimes(2))
fprintf('Comparisons count: %i\n',length(structComparisons))
fprintf('Regions count: %i\n',sum(str2double(switchRegions(:,2))))
fprintf('Permutation cycles: %i\n',intCycles)
fprintf('Permutation alpha: %.2f\n',dblAlpha)
fprintf('T-test type: %s ',strTestType)
if boolUseSingleTail
    fprintf('1-tailed\n')
else
    fprintf('2-tailed\n')
end
fprintf('\n')
if boolUseSubjectExclusion
    fprintf('Warning: subject exclusion is on\n')
end
if ~exist(dirOutput,'dir')
    fprintf('Warning: output folder does not exist and will be automatically created\n')
end
fprintf('\n')

%confirm user-input
strUserInput = input('Continue with the above settings? (Y/N): ','s');
if isempty(strfind(strUserInput,'Y')) && isempty(strfind(strUserInput,'y'))
    error('Operation aborted')
end

%validate output directory
if ~exist(dirOutput,'dir')
    [boolSuccess,strMessage,strMessageID] = mkdir(dirOutput);
    if ~boolSuccess
        error('%s: %s',strMessageID,strMessage)
    end
end

%load in all processed wavelet data MAT files
fprintf('Loading inputs:\n')
for i = 1:length(lstInfiles)
    fprintf('\t%s\n',lstInfiles{i})
    load([dirInfolder,lstInfiles{i}])
    
    %save into data structure
    dataWavelet(i) = Data;
    if exist('History','var')
        dataHistory(i) = History;
    end
    
    clear Data History
end
fprintf('\n')
clear i lstInfiles

%regions definitions-------------------------------------------------------

load(pathElectrodeMap)

structRegions(1).name = 'all-PFC';
structRegions(1).electrodes = cat(2, wholehead_elecs{1,2},wholehead_elecs{2,2}, wholehead_elecs{3,2});
structRegions(2).name = 'left-PFC';
structRegions(2).electrodes = wholehead_elecs{2,2};
structRegions(3).name = 'mid-PFC';
structRegions(3).electrodes = wholehead_elecs{1,2};
structRegions(4).name = 'right-PFC';
structRegions(4).electrodes = wholehead_elecs{3,2};
structRegions(5).name = 'all-Central';
structRegions(5).electrodes = cat(2, wholehead_elecs{8,2},wholehead_elecs{4,2}, wholehead_elecs{9,2});
structRegions(6).name = 'left-Central';
structRegions(6).electrodes = wholehead_elecs{8,2};
structRegions(7).name = 'mid-Central';
structRegions(7).electrodes = wholehead_elecs{4,2};
structRegions(8).name = 'right-Central';
structRegions(8).electrodes = wholehead_elecs{9,2};
structRegions(9).name = 'all-Post';
structRegions(9).electrodes = cat(2, wholehead_elecs{6,2},wholehead_elecs{5,2}, wholehead_elecs{7,2});
structRegions(10).name = 'left-Post';
structRegions(10).electrodes = wholehead_elecs{6,2};
structRegions(11).name = 'mid-Post';
structRegions(11).electrodes = wholehead_elecs{5,2};
structRegions(12).name = 'right-Post';
structRegions(12).electrodes = wholehead_elecs{7,2};

clear wholehead_elecs

%begin processing algorithm------------------------------------------------

fprintf('\n')

%loop through all comparisons
for structComparison = structComparisons
    
    fprintf('Loading data for %s %s Vs %s %s ...', ...
    	structComparison.pos.dose,structComparison.pos.cond,structComparison.neg.dose,structComparison.neg.cond)
    
    %find the data in the processed wavelet data
    idxDosage = find(cellfun(@(x) strcmp(x,structComparison.pos.dose),{dataWavelet.field}));
    idxCondition = find(cellfun(@(x) strcmp(x,structComparison.pos.cond),{dataWavelet.subfield}));
    idxGroup1 = intersect(idxDosage,idxCondition);
    idxDosage = find(cellfun(@(x) strcmp(x,structComparison.neg.dose),{dataWavelet.field}));
    idxCondition = find(cellfun(@(x) strcmp(x,structComparison.neg.cond),{dataWavelet.subfield}));
    idxGroup2 = intersect(idxDosage,idxCondition);
    clearvars idxDosage idxCondition

    %sanity check: make sure subject numbers match
    if ~isequal(dataWavelet(idxGroup1).subject_list,dataWavelet(idxGroup2).subject_list)
        error('subject entries do not match, data is not paired: %s %s Vs %s %s', ...
            structComparison.pos.dose,structComparison.pos.cond,structComparison.neg.dose,structComparison.neg.cond)
    end
    
    %get the data for each group <subject X electrode X frequency X time>
    dataGroup1 = dataWavelet(idxGroup1).matrix;
    dataGroup2 = dataWavelet(idxGroup2).matrix;
    
    %sanity check: make sure both data matrices have the same dimension
    if ~isequal(size(dataGroup1),size(dataGroup2))
        error('data matrices do not match, data is not paired: %s %s Vs %s %s', ...
            structComparison.pos.dose,structComparison.pos.cond,structComparison.neg.dose,structComparison.neg.cond)
    end

    %check if subjection exclusion is enabled
    if boolUseSubjectExclusion
        %find the index of subjects and remove them
        idxDelete = find(ismember(dataWavelet(idxGroup1).subject_list,lstExcludedSubjects));
        dataGroup1(idxDelete,:,:,:) = [];
        idxDelete = find(ismember(dataWavelet(idxGroup2).subject_list,lstExcludedSubjects));
        dataGroup2(idxDelete,:,:,:) = [];
        
        clearvars idxDelete
    end
    
    fprintf('done\n\n')
    
    %loop through all regions(keep)
    for cellRegion = switchRegions' %transpose needed to iterate by rows
        %check if this region should be ran
        if ~str2double(cellRegion{2})
            fprintf('Region:%14s | skipped\n\n',cellRegion{1})
            continue
        end
        
        %find the region to run
        idxRegion = find(cellfun(@(x) strcmp(x,cellRegion{1}),{structRegions.name}));
        
        %extract data from region of interest
        dataGroup1Subset = dataGroup1(:,structRegions(idxRegion).electrodes,:,:);
        dataGroup2Subset = dataGroup2(:,structRegions(idxRegion).electrodes,:,:);
        clearvars idxRegion
                
        %average down the region to 1 giant electrode <subject X frequency X time>
        dataGroup1Subset = squeeze(nanmean(dataGroup1Subset,2));
        dataGroup2Subset = squeeze(nanmean(dataGroup2Subset,2));
        
        %iterate through all frequency and time points
        dataResult.all_t = zeros(size(dataGroup1Subset,2),size(dataGroup1Subset,3));	%used to be called all_t_matrix, in t-values
        dataResult.all_p = zeros(size(dataGroup1Subset,2),size(dataGroup1Subset,3));	%new, in p-values
        dataResult.sig_t = zeros(size(dataGroup1Subset,2),size(dataGroup1Subset,3));	%used to be called sig_matrix, in t-values
        dataResult.sig_p = zeros(size(dataGroup1Subset,2),size(dataGroup1Subset,3));	%new, in p-values
        for i = 1:size(dataGroup1Subset,2)  %frequency
            fprintf('Region:%14s | Frequency:%3i | Elapsed: ',cellRegion{1},arrFrequency(1) + i - 1)
            tic
            
            for j = 1:size(dataGroup1Subset,3)    %time
                %extract subject data
                arrGroup1 = dataGroup1Subset(:,i,j);
                arrGroup2 = dataGroup2Subset(:,i,j);
                
                %call legacy function
                [dblPValue,dblTValue] = eval([strFunction,'(arrGroup1,arrGroup2,strTestType,intCycles)']);
                
                clearvars arrGroup1 arrGroup2
                
                %save the result in all data matrix
                dataResult.all_t(i,j) = dblTValue;
                dataResult.all_p(i,j) = dblPValue;
                
                %perform significance comparison
                if boolUseSingleTail
                    %assume right tail only
                    if (dblPValue > (1-dblAlpha))
                        dataResult.sig_t(i,j) = dblTValue;
                        dataResult.sig_p(i,j) = dblPValue;
                    end
                else
                    if (dblPValue < dblAlpha/2) || (dblPValue > (1-dblAlpha/2))
                        dataResult.sig_t(i,j) = dblTValue;
                        dataResult.sig_p(i,j) = dblPValue;
                    end
                end
                
                clearvars dblPValue dblTValue
            end 
            
            fprintf('%.3f seconds\n',toc)
        end
        
        fprintf('\n')
        
        %create generic variable to store results
        Data = dataResult;
        Data.phase = dataWavelet(idxGroup1).phase;
        Data.type = dataWavelet(idxGroup1).type;
        Data.frequency = dataWavelet(idxGroup1).frequency;
        Data.time = dataWavelet(idxGroup1).time;
        Data.label = structComparison.name;
        Data.region = cellRegion{1};
        clearvars dataResult
        
        %save history information
        if exist('dataHistory','var')
            myHistory.source.Group1_name = dataHistory(idxGroup1).source_name;
            myHistory.source.Group1_id = dataHistory(idxGroup1).source_id;
            myHistory.source.Group2_name = dataHistory(idxGroup2).source_name;
            myHistory.source.Group2_id = dataHistory(idxGroup2).source_id;
        end
        History = myHistory;
        
        %create save filename
        strBuilder = [strFilePrefix,'_'];
        strBuilder = [strBuilder,structComparison.name,'_'];    %comparison name
        strBuilder = [strBuilder,cellRegion{1}];                %region name
        
        %save the file
        fprintf('Saving...')
        save([dirOutput,strBuilder],'Data','History')
        fprintf('done\n\n')
        
        clearvars Data strBuilder dataGroup1Subset dataGroup2Subset
    end
    
    clearvars dataGroup1 dataGroup2
end

clearvars

fprintf('Done\n\n')
beep

%% Plot
%
%This makes plots for all data (all matrices and significant matrices),
%mainly used for debugging. For actual data plotting use multiplot below.
%Note that for sig p values, the 0's are adjusted to 0.5 to make it colored
%the same as sig t values

%Note: load in result MAT file of interest

%user-defined vars
arrColorScale_t = [-4 4];   %for t-values
arrColorScale_p = [0 0.1];    %for p-values
boolAdjustPValues = true;   %this makes all p-values that passed the p > (1-alpha) test to be flipped about 0.5 to be within the region of [0 alpha]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%begin algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%data validation
if ~exist('Data','var') || ~exist('History','var')
    error('Input data not found')
end

%get time and frequency arrays
arrFrequency = Data.frequency;
arrTime = Data.time;

%p-value adjustment
if boolAdjustPValues
    for i = 1:size(Data.all_p,1)    %loop through all frequencies
        %find indexes to fix
        idxRightTails_all = find(Data.all_p(i,:) > 0.5);
        idxRightTails_sig = find(Data.sig_p(i,:) > 0.5);
        idxRightTails_sig2 = find(Data.sig_p(i,:) == 0);
        
        %adjust the values
        Data.all_p(i,idxRightTails_all) = 1 - Data.all_p(i,idxRightTails_all);
        Data.sig_p(i,idxRightTails_sig) = 1 - Data.sig_p(i,idxRightTails_sig);
        Data.sig_p(i,idxRightTails_sig2) = 0.05;
    end
end

%create plot
figure
subplot(2,2,1), ...
    tftopo(Data.all_t,arrTime,arrFrequency,'title','all t','verbose','off'), ...
    caxis(arrColorScale_t), ...
    colorbar, ...
subplot(2,2,2), ...
    tftopo(Data.sig_t,arrTime,arrFrequency,'title','sig t','verbose','off'), ...
    caxis(arrColorScale_t), ...
    colorbar, ...
subplot(2,2,3), ...
    tftopo(Data.all_p,arrTime,arrFrequency,'title','all p','verbose','off'), ...
    caxis(arrColorScale_p), ...
    colorbar, ...
subplot(2,2,4), ...
    tftopo(Data.sig_p,arrTime,arrFrequency,'title','sig p','verbose','off'), ...
    caxis(arrColorScale_p), ...
    colorbar;
% ax = axes('Units','Normal','Position',[0.3 -0.85 .85 .85],'Visible','off');
% set(get(ax,'Title'),'Visible','on');
% title(History.vars.dirInfolder)
%% Multiplot
%
%This makes all plots in the outfolder under the current directory.
%
%Note: run this code cell after setting the working path to the root of the
%permutation files
%
%Changelog:
%   3/28/2013:  Added processing of History variable if it exists. If
%               multiple files exist for the same data type (i.e. if
%               subject groups were processed in different batches), then
%               the History from the first file will be used.

clear,clc

%user-defined vars
arrColorScale_t = [-4 4];   %for t-values
arrColorScale_p = [0 0.1];    %for p-values
boolAdjustPValues = true;   %this makes all p-values that passed the p > (1-alpha) test to be flipped about 0.5 to be within the region of [0 alpha]
boolPublicationMode = false; %this disables titles and source labeling in the graphs

%task-specific vars
strPublicationFolder = 'Publication/';	%this is the subfolder that would be created within the normal output directory in publication mode

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%begin algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%find all files to process
lstFiles = dir('PRMT*');
lstFiles = {lstFiles.name};

%check if there are files to process
if isempty(lstFiles)
    error('Cannot find files to process with prefix: %s',strFilePrefix)
end

%create output folders
if boolPublicationMode
    dirOutput_all_t = [strPublicationFolder,'All t figures/'];
    dirOutput_sig_t = [strPublicationFolder,'Significant t figures/'];
    dirOutput_all_p = [strPublicationFolder,'All p figures/'];
    dirOutput_sig_p = [strPublicationFolder,'Significant p figures/'];
    dirOutput_combined = '';    %combined figures for investigational purposes only
    
    if ~exist(dirOutput_all_t,'dir')
        mkdir(dirOutput_all_t);
    end
    if ~exist(dirOutput_sig_t,'dir')
        mkdir(dirOutput_sig_t);
    end
    if ~exist(dirOutput_all_p,'dir')
        mkdir(dirOutput_all_p);
    end
    if ~exist(dirOutput_sig_p,'dir')
        mkdir(dirOutput_sig_p);
    end
else
    dirOutput_all_t = 'All t figures/';
    dirOutput_sig_t = 'Significant t figures/';
    dirOutput_all_p = 'All p figures/';
    dirOutput_sig_p = 'Significant p figures/';
    dirOutput_combined = 'Combined figures/';
    
    if ~exist(dirOutput_all_t,'dir')
        mkdir(dirOutput_all_t);
    end
    if ~exist(dirOutput_sig_t,'dir')
        mkdir(dirOutput_sig_t);
    end
    if ~exist(dirOutput_all_p,'dir')
        mkdir(dirOutput_all_p);
    end
    if ~exist(dirOutput_sig_p,'dir')
        mkdir(dirOutput_sig_p);
    end
    if ~exist(dirOutput_combined,'dir')
        mkdir(dirOutput_combined);
    end
end

%process all files
figure
for cellFile = lstFiles
    strFullFilename = cell2mat(cellFile);
    [~,strFilename,~] = fileparts(strFullFilename);
    
    load(strFullFilename)
    
    %data validation
    if ~exist('Data','var') || ~exist('History','var')
        error('Input data not found: %s',strFullFilename)
    end

    %p-value adjustment
    if boolAdjustPValues
        for i = 1:size(Data.all_p,1)    %loop through all frequencies
            %find indexes to fix
            idxRightTails_all = find(Data.all_p(i,:) > 0.5);
            idxRightTails_sig = find(Data.sig_p(i,:) > 0.5);
            idxRightTails_sig2 = find(Data.sig_p(i,:) == 0);

            %adjust the values
            Data.all_p(i,idxRightTails_all) = 1 - Data.all_p(i,idxRightTails_all);
            Data.sig_p(i,idxRightTails_sig) = 1 - Data.sig_p(i,idxRightTails_sig);
            Data.sig_p(i,idxRightTails_sig2) = 0.05;
        end
    end
    
    %build title
    strBuilder = Data.phase;
    strBuilder = [strBuilder,'-',Data.type];
    strBuilder = [strBuilder,'-',num2str(min(Data.frequency)),'-',num2str(max(Data.frequency)),'Hz'];
    strBuilder = [strBuilder,'-',strrep(Data.label,'_','-')];
    strBuilder = [strBuilder,'-',Data.region];

    %build source string
    if ismember('source',fieldnames(History))
        strSource1_name = History.source.Group1_name;
        strSource1_id = History.source.Group1_id;
        strSource2_name = History.source.Group1_name;
        strSource2_id = History.source.Group2_id;
        if (strcmp(strSource1_name,strSource2_name) && strcmp(strSource1_id,strSource2_id))
            strSource = ['Source: ',strrep(strSource1_name,'_','\_'),' (',strSource1_id,')'];
        else
            strSource = {['Source: ',strrep(strSource1_name,'_','\_'),' (',strSource1_id,')'],[strrep(strSource2_name,'_','\_'),' (',strSource2_id,')']};
        end
    else
        strSource = ['Source: ',History.vars.dirInfolder];
    end
    
    %plot individual figures
    
    %all t
    if boolPublicationMode
        strTitle = ' ';
    else
        strTitle = [strBuilder,' {\bf{\color{blue}All-t}}'];
    end
    h = tftopo(Data.all_t,Data.time,Data.frequency,'title',strTitle,'verbose','off');
    caxis(arrColorScale_t), colorbar;
    if ~boolPublicationMode
        lblText = axes('Units','Normal','Position',[0.35 -0.86 .85 .85],'Visible','off');
        set(get(lblText,'Title'),'Visible','on'), title(strSource,'FontSize',6);
    end
    saveas(gcf,[dirOutput_all_t,strFilename,'.tif'],'tiff');
    clf
    
    %sig t
    if boolPublicationMode
        strTitle = ' ';
    else
        strTitle = [strBuilder,' {\bf{\color{blue}Sig-t}}'];
    end
    h = tftopo(Data.sig_t,Data.time,Data.frequency,'title',strTitle,'verbose','off');
    caxis(arrColorScale_t), colorbar;
    if ~boolPublicationMode
        lblText = axes('Units','Normal','Position',[0.35 -0.86 .85 .85],'Visible','off');
        set(get(lblText,'Title'),'Visible','on'), title(strSource,'FontSize',6);
    end
    saveas(gcf,[dirOutput_sig_t,strFilename,'.tif'],'tiff');
    clf
    
    %all p
    if boolPublicationMode
        strTitle = ' ';
    else
        strTitle = [strBuilder,' {\bf{\color{blue}All-p}}'];
    end
    h = tftopo(Data.all_p,Data.time,Data.frequency,'title',strTitle,'verbose','off');
    caxis(arrColorScale_p), colorbar;
    if ~boolPublicationMode
        lblText = axes('Units','Normal','Position',[0.35 -0.86 .85 .85],'Visible','off');
        set(get(lblText,'Title'),'Visible','on'), title(strSource,'FontSize',6);
    end
    saveas(gcf,[dirOutput_all_p,strFilename,'.tif'],'tiff');
    clf
    
    %sig p
    if boolPublicationMode
        strTitle = ' ';
    else
        strTitle = [strBuilder,' {\bf{\color{blue}Sig-p}}'];
    end
    h = tftopo(Data.sig_p,Data.time,Data.frequency,'title',strTitle,'verbose','off');
    caxis(arrColorScale_p), colorbar;
    if ~boolPublicationMode
        lblText = axes('Units','Normal','Position',[0.35 -0.86 .85 .85],'Visible','off');
        set(get(lblText,'Title'),'Visible','on'), title(strSource,'FontSize',6);
    end
    saveas(gcf,[dirOutput_sig_p,strFilename,'.tif'],'tiff');
    clf
    
    %combined
    if ~boolPublicationMode
        subplot(2,2,1), ...
            h = tftopo(Data.all_t,Data.time,Data.frequency,'title','All t','verbose','off');
            caxis(arrColorScale_t), colorbar;
        subplot(2,2,2), ...
            h = tftopo(Data.sig_t,Data.time,Data.frequency,'title','Sig t','verbose','off');
            caxis(arrColorScale_t), colorbar;
        subplot(2,2,3), ...
            h = tftopo(Data.all_p,Data.time,Data.frequency,'title','All p','verbose','off');
            caxis(arrColorScale_p), colorbar;
        subplot(2,2,4), ...
            h = tftopo(Data.sig_p,Data.time,Data.frequency,'title','Sig p','verbose','off');
            caxis(arrColorScale_p), colorbar;
            lblText = axes('Units','Normal','Position',[0.35 -0.86 .85 .85],'Visible','off');
            set(get(lblText,'Title'),'Visible','on'), title(strSource,'FontSize',6);
        saveas(gcf,[dirOutput_combined,strFilename,'.tif'],'tiff');
        clf
    end
    
    clear Data History
end

close,beep