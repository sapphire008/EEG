%Extract data for statistical analysis
%Author: Julian Y. Cheng
%
%Note: LOAD IN DESIRED ERSP/ITC FIRST (do not load both at the same time)

%1/4/2013 run
%Frequency bands: theta(4-8),alpha(8-15),beta(15-30)
%Time interval: 0-200ms
%Conditions: Plac contrast, Drug contrast
%Regions: left PFC, mid PFC, right PFC, all PFC

%1/10/2013 run
%Mode: MAX
%Frequency bands: theta(4-8),alpha(8-15),beta(15-30)
%Time interval: 0-250ms
%Conditions: Plac contrast, Drug contrast, Qmoney contrast
%Regions: left PFC, mid PFC, right PFC, all PFC
%Notes: This is the first implementation of flagMode, and an exception is
%       created to construct "Qmoney contrast" on-the-fly; "Qmoney green"
%       should work also, but needs to be coded as "Qmoney cong" for the 
%       algorithm to recognize the variables

%% Extract data from ERSP/ITC

%user-specific vars
flagMode = 1;   %0 = find average, 1 = find latency to (first) max value
lstFrequencies = [{[4 8]},{[8 15]},{[15 30]}];                              %each cell is a different frequency pair
lstTimeIntervals = [{[0 250]}];                                             %same format as lstFrequencies
lstConditions = [{[{'Plac'},{'contrast'}]}, ...                             %each cell has a cell arrays of condition pairs
                {[{'Drug'},{'contrast'}]}, ...
                {[{'Qmoney'},{'contrast'}]}];
lstRegions = [{'left_PFC'},{'mid_PFC'},{'right_PFC'},{'all_PFC'}];          %same format as lstFrequencies
strFilename = 'theta_alpha_beta_0-250ms_PFC';             %will append mode specifier string to beginning of filename

%task-specific vars
strPrefix0 = 'AVG'; %for mode 0
strPrefix1 = 'MAX'; %for mode 1
dirOutput = '/nfs/erp-modaf/mc/MC-epoched-Feb11-Glenn/Cue/Statistics/';

%task-specific vars (electrode locations)
load('/nfs/erp-modaf/elec_files/wholehead_elecs.mat');
left_PFC = wholehead_elecs{2,2};
mid_PFC = wholehead_elecs{1,2};
right_PFC = wholehead_elecs{3,2};
left_Central = wholehead_elecs{8,2};
mid_Central = wholehead_elecs{4,2};
right_Central = wholehead_elecs{9,2};
left_Post = wholehead_elecs{6,2};
mid_Post = wholehead_elecs{5,2};
right_Post = wholehead_elecs{7,2};
all_PFC = cat(2, wholehead_elecs{1,2},wholehead_elecs{2,2}, wholehead_elecs{3,2});
all_Central = cat(2, wholehead_elecs{8,2},wholehead_elecs{4,2}, wholehead_elecs{9,2});
all_Post = cat(2, wholehead_elecs{6,2},wholehead_elecs{5,2}, wholehead_elecs{7,2});
all_Head = 1:128;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%begin algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc

fprintf('Data extraction for statistical analysis\nAuthor: Julian Y. Cheng\n\n');
fprintf('Processing mode is %s\n\n',eval(['strPrefix',num2str(flagMode)]))

%input data integrity check
if ~exist('ERSP','var') && ~exist('ITC','var')
    error('Failed to find variable: ERSP or ITC')
elseif exist('ERSP','var') && exist('ITC','var')
    error('Cannot have ERSP and ITC concurrently in the workspace')
elseif ~exist('freqs','var')
    error('Failed to find variable: freqs')
elseif ~exist('timesout','var')
    error('Failed to find variable: timesout')
end
if exist('ERSP','var')
    strDataType = 'ERSP';
    myData = ERSP;
    boolIsERSP = true;
else
    strDataType = 'ITC';
    myData = ITC;
    boolIsERSP = false;
end

%construct history
arrTime = fix(clock);
%History.script = mfilename('fullpath'); %doesn't work in cell mode!
History.mode = eval(['strPrefix',num2str(flagMode)]);
History.date = [num2str(arrTime(2)),'/',num2str(arrTime(3)),'/',num2str(arrTime(1)),' ',num2str(arrTime(4)),':',num2str(arrTime(5))];
History.vars.lstFrequencies = lstFrequencies;
History.vars.lstTimeIntervals = lstTimeIntervals;
History.vars.lstConditions = lstConditions;
History.vars.lstRegions = lstRegions;
History.vars.dirOutput = dirOutput;
History.vars.strDataType = strDataType;

%process all combinations of input variables
structData = struct();  %holds all data combinations
for cellFrequencyRange = lstFrequencies
    arrFrequencyRange = cell2mat(cellFrequencyRange);
    
    for cellTimeInterval = lstTimeIntervals
        arrTimeInterval = cell2mat(cellTimeInterval);
        
        for cellCondition = lstConditions
            strDosage = cellCondition{1}{1};
            strCondition = cellCondition{1}{2};
            
            for cellRegion = lstRegions
                strRegion = cell2mat(cellRegion);
                
                fprintf('Processing %2i-%2iHz %i-%ims %s_%s %s\t...',arrFrequencyRange(1),arrFrequencyRange(2),arrTimeInterval(1),arrTimeInterval(2),strDosage,strCondition,strRegion)
                
                %find frequency and time interval indexes
                arrFrequencyIndexes = find(freqs >= arrFrequencyRange(1),1):find(freqs <= arrFrequencyRange(2),1,'last');
                arrTimeIntervalIndexes = find(timesout >= arrTimeInterval(1),1):find(timesout <= arrTimeInterval(2),1,'last');
                
                %set main processing loop condition
                %will execute the processing loop twice for Qmoney
                if strcmp(strDosage,'Qmoney')   
                    %check if required fields exist in myData
                    if ~isfield(myData,'Drug') || ~isfield(myData,'Plac')
                        error('Qmoney processing requires both Drug and Plac fields to be present in ERSP')
                    end
                        
                    boolIsQmoney = true;
                else
                    boolIsQmoney = false;
                end
                
                for iterator = 1:2  %this is to process both conditions for Qmoney
                    %get the data of interest
                    if boolIsQmoney   %exception: prepare to construct Qmoney-related conditions on-the-fly
                        if (iterator == 1)  %process drug first, then plac
                            structAllSubjectData = myData.Drug;
                        else
                            structAllSubjectData = myData.Plac;
                        end
                    else
                        eval(['structAllSubjectData = myData.',strDosage,';']);
                    end

                    %loop through all subjects to extract frequency X time data (also saves which subjects are in the dataset)
                    lstSubjects = [];
                    matrixAllSubjectData = [];  %saves data for all subjects across electrodes, frequencies, time
                    for idxSubject = 1:length(structAllSubjectData)
                        %get the current subject data
                        structSubject = structAllSubjectData{idxSubject};

                        %check if this is empty
                        if isempty(structSubject)
                            continue
                        else
                            lstSubjects(end+1) = idxSubject;
                        end

                        %get the electrode data of interest
                        eval(['arrAllElectrodeData = structSubject.',strCondition,'.lead;'])

                        %process electrodes that are within region of interest
                        for intElectrode = all_Head
                            %check if electrode exists
                            if (intElectrode > length(arrAllElectrodeData))
                                matrixAllSubjectData(length(lstSubjects),intElectrode,1:length(arrFrequencyIndexes),1:length(arrTimeIntervalIndexes)) = NaN;
                                continue
                            end

                            %get current electrode data
                            matrixElectrode = arrAllElectrodeData{intElectrode};

                            %exclude if:
                            %1. electrode is not in region of interest
                            %2. electrode has no data
                            %3. electrode data has an average of 0
                            if eval(['~ismember(intElectrode,',strRegion,')']) || ...
                               isempty(matrixElectrode) || ...
                               (mean2(matrixElectrode) == 0)
                                matrixAllSubjectData(length(lstSubjects),intElectrode,1:length(arrFrequencyIndexes),1:length(arrTimeIntervalIndexes)) = NaN;
                                continue
                            end

                            %data is valid, save
                            if boolIsERSP
                                matrixAllSubjectData(length(lstSubjects),intElectrode,1:length(arrFrequencyIndexes),1:length(arrTimeIntervalIndexes)) = ...
                                    matrixElectrode(arrFrequencyIndexes,arrTimeIntervalIndexes);
                            else
                                matrixAllSubjectData(length(lstSubjects),intElectrode,1:length(arrFrequencyIndexes),1:length(arrTimeIntervalIndexes)) = ...
                                    abs(matrixElectrode(arrFrequencyIndexes,arrTimeIntervalIndexes));
                            end
                        end
                    end
                    
                    if ~boolIsQmoney
                        %this is not Qmoney, exit now
                        break
                    elseif (iterator == 1)
                        %this is drug iteration for Qmoney, save the result
                        matrixAllSubjectData_Drug = matrixAllSubjectData;
                    else
                        %this is plac iteration for Qmoney, construct Qmoney matrix
                        matrixAllSubjectData_Plac = matrixAllSubjectData;
                        matrixAllSubjectData = matrixAllSubjectData_Drug - matrixAllSubjectData_Plac;
                    end
                end
                
                switch flagMode
                    case 0
                        %average across electrode,frequency,time
                        %note that the original dimensions are <subject X electrode X frequency X time>
                        matrixAllSubjectData = squeeze(nanmean(matrixAllSubjectData,2));    %across electrode (yeilds 1 giant electrode)
                        matrixAllSubjectData = squeeze(nanmean(matrixAllSubjectData,2));    %across frequency (compresses the frequency band)
                        arrAllSubjectData = squeeze(nanmean(matrixAllSubjectData,2));       %across time
                    case 1
                        %average across electrode
                        matrixAllSubjectData = squeeze(nanmean(matrixAllSubjectData,2));    %across electrode (yeilds 1 giant electrode)
                        
                        %find maximum for each subject
                        arrAllSubjectData = zeros(size(matrixAllSubjectData,1),1);
                        for idxSubject = 1:size(matrixAllSubjectData,1)
                            matrixThisSubject = squeeze(matrixAllSubjectData(idxSubject,:,:));
                            [~,idxMax] = max(matrixThisSubject(:));                     %note: this is a linear index, as if the entire matrix was transformed into 1 giant array
                            [idxRow,idxCol] = ind2sub(size(matrixThisSubject),idxMax);  %row index = frequency (ignored), col index = time
                            arrAllSubjectData(idxSubject) = timesout(arrTimeIntervalIndexes(idxCol));    %this is to reconstruct the time value because the index we found is from an extracted set
                        end
                    otherwise
                        error('Invalid flag mode specified')
                end
                
                %save the data
                if isempty(fieldnames(structData))
                    structData(1).type.frequency = arrFrequencyRange;
                    structData(1).type.time = arrTimeInterval;
                    structData(1).type.dose = strDosage;
                    structData(1).type.condition = strCondition;
                    structData(1).type.region = strRegion;
                    structData(1).subjects = lstSubjects;
                    structData(1).array = arrAllSubjectData;
                else
                    structData(end+1).type.frequency = arrFrequencyRange;
                    structData(end).type.time = arrTimeInterval;
                    structData(end).type.dose = strDosage;
                    structData(end).type.condition = strCondition;
                    structData(end).type.region = strRegion;
                    structData(end).subjects = lstSubjects;
                    structData(end).array = arrAllSubjectData;
                end
                
                fprintf('done\n')
            end
        end
    end
end

%save the variables
fprintf('\nSaving...')
save([dirOutput,eval(['strPrefix',num2str(flagMode)]),'_',strFilename],'structData','History')
fprintf('done\n\n')

%% Export into csv

%Note: this outputs 1 csv file for every data array, so combine them in
%Windows however you want later

%task-specific vars
dirOutput = '/nfs/erp-modaf/mc/MC-epoched-Feb11-Glenn/Cue/Statistics/Output/';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%begin algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check for input data
if ~exist('structData','var')
    error('Failed to find input data')
end

%process all data arrays
for idxSubject = 1:length(structData)
    
    %build save filename
    myInfo = structData(idxSubject).type;
    strPrefix = History.mode;
    strBuilder = [strPrefix,'_'];
    strBuilder = [strBuilder,num2str(myInfo.frequency(1)),'-',num2str(myInfo.frequency(2)),'Hz_'];
    strBuilder = [strBuilder,num2str(myInfo.time(1)),'-',num2str(myInfo.time(2)),'ms_'];
    strBuilder = [strBuilder,myInfo.dose,'_',myInfo.condition,'_',myInfo.region,'.csv'];
    
    %build data matrix <subject id X data>
    matrixOutput(:,1) = structData(idxSubject).subjects';
    matrixOutput(:,2) = structData(idxSubject).array;
    
    %save file
    csvwrite([dirOutput,strBuilder],matrixOutput)
    
    clear myInfo strPrefix strBuilder matrixOutput
end