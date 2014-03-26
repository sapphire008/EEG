%% bootstrap ERSP (adapted from Josh's version)
%This file is a modified version of CMET_bootstrap_corrected_between_gpu.m
%from the SC study, which was originally given to Julian Cheng by Josh
%Phillips
%
%Changelog:
%12/19/2012: modified for use with MC/PC headplots

%LOAD INPUT 4D HDPLT MAT FILES FROM srPlotHeadplots.m FIRST
%MAKE SURE THIS IS RUN ON A MACHINE WITH GPU HARDWARE

%% Perform permutations

addpath /home/phillips/Scripts/
addpath /nfs/pkg64/gpu/

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
fcz = [19];
all_PFC = cat(2, wholehead_elecs{1,2},wholehead_elecs{2,2}, wholehead_elecs{3,2});
all_Central = cat(2, wholehead_elecs{8,2},wholehead_elecs{4,2}, wholehead_elecs{9,2});
all_Post = cat(2, wholehead_elecs{6,2},wholehead_elecs{5,2}, wholehead_elecs{7,2});
all_Head = 1:128;

%user-defined vars

%select comparison group to run
% 0 = Drug Red vs Green
% 1 = Plac Red vs Green
% 2 = Q money, i.e. Drug contrast vs Plac contrast
% 3 = Plac Red vs Null (all zeros)
% (anything else) = Other, makes use of lstCustomGroup
flagGroup = 1;
lstCustomGroup = [{'Drug_incong'} {'Drug_cong'}];
strGroup = 'Drug_incongVcong';  %used in filename

%select regions to run
region_average_names = [{'all_Head'}];

%permutation parameters
perm_number = 50000;
alpha = .05;
intRejectLimit = 17;    %rejects electrode if excluded subjects exceed this number
    
%outfile saving
dirOutput = '/nfs/erp-modaf/mc/MC-epoched-Feb11-Glenn/Cue/Headplots/Bootstrap/';
boolDoNotOverwrite = 1; %will warn and append numbers to new filename

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%begin algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%auto path generation, used if input mat file contains variable strPhase
if exist('strPhase','var')
    dirOutput = strrep(dirOutput,'Cue',strPhase);
end

switch flagGroup
    case 0
        group = [{'Drug_incong'} {'Drug_cong'}]; %Drug Red vs Green
        strGroup = 'Drug_RvG';
    case 1
        group = [{'Plac_incong'} {'Plac_cong'}]; %Plac Red vs Green
        strGroup = 'Plac_RvG';
    case 2
        group = [{'Drug_contrast'} {'Plac_contrast'}]; %Q money
        strGroup = 'Qmoney';
    case 3
        group = [{'Plac_cong'} {'Null'}]; %QPlac Red vs Null (all zeros)
        matrixSubjectByElectrode_Null = zeros(size(matrixSubjectByElectrode_Plac_cong));
        matrixSubjectByElectrode_Null(isnan(matrixSubjectByElectrode_Plac_cong)) = NaN;
        strGroup = 'Plac_GvN';
    otherwise
        group = lstCustomGroup;
end

%get input data matrices
eval([group{1},' = matrixSubjectByElectrode_',group{1},';'])
eval([group{2},' = matrixSubjectByElectrode_',group{2},';'])

for l = 1:length(region_average_names)
    
    eval(['region_average_elecs =', region_average_names{l},';'])
    region_average_elecs = sort(region_average_elecs);

    elecs = 1;

    % squeeze and nanmean across region to create a single electrode
    % (average)
    eval(['first_group_contrast_region = ',group{1},'(:,region_average_elecs);'])
    eval(['second_group_contrast_region = ',group{2},'(:,region_average_elecs);'])

    intMaxNumberElectrodes = length(region_average_elecs);
    intElectrodes = size(eval([group{1}]),2);
    
    actual_t_matrix = zeros(1,intElectrodes);
    sig_matrix = zeros(1,intElectrodes);
    lstPValues = NaN(1,intElectrodes);
    lstExcludedElectrodes = zeros(1,intElectrodes);

    tic
    for idxElectrode = 1:intMaxNumberElectrodes

        cond1_vect = first_group_contrast_region(:,idxElectrode);
        cond2_vect = second_group_contrast_region(:,idxElectrode);         

        combined_vector = [cond1_vect ; cond2_vect];

        %remove NaNs
        idxToRemove = isnan(combined_vector);
        combined_vector(idxToRemove) = [];

        fraction = bootstrap(combined_vector, length(cond1_vect), length(cond2_vect), perm_number);
        
        clear combined_vector group_* sequence

        sample_mean_difference = nanmean(cond1_vect)-nanmean(cond2_vect);

        idxActualElectrodePos = find(ismember(1:intElectrodes,region_average_elecs(idxElectrode)));

        if (length(find(idxToRemove))/2 <= intRejectLimit)
            actual_t_matrix(idxActualElectrodePos) = sample_mean_difference;

            if fraction < alpha/2 || fraction > (1 - alpha/2)
                sig_matrix(idxActualElectrodePos) = sample_mean_difference;
            else
                sig_matrix(idxActualElectrodePos)= 0;
            end
            
            strElectrodeRejected = '';
        else
            actual_t_matrix(idxActualElectrodePos) = NaN;
            sig_matrix(idxActualElectrodePos)= NaN;
            
            strElectrodeRejected = 'REJECTED';
        end
        
        lstPValues(idxElectrode) = fraction;
        lstExcludedElectrodes(idxElectrode) = length(find(idxToRemove))/2;

        clear h p ci stats upper_t lower_t pseudo* cond1_vect cond2_vect combined_vector sorted*

        disp(['Region: ', region_average_names{l},' | Elec: ',num2str(idxElectrode),'/',num2str(intMaxNumberElectrodes),' | Missing ',num2str(length(find(idxToRemove))/2),' subject(s)  ',strElectrodeRejected])
    end
    toc
    
    %build save filename
    strFilePostfix = strFilename(strfind(strFilename,strPhase)+length(strPhase)+1:end-4);
    pathOutfile = [dirOutput,'HDPLT_BSTRAP_',strFilePostfix,'_',strGroup,'_',region_average_names{l},'.mat'];
    
    strSourceFile = [dirOutput,strFilename];
    
    %non-overwrite protection
    intCounter = 2;
    while boolDoNotOverwrite && exist(pathOutfile,'file')
        if isempty(regexp(pathOutfile,'_\d+.mat','match'))
            %add number to end
            pathOutfile = strrep(pathOutfile,'.mat',['_',num2str(intCounter),'.mat']);
            intCounter = intCounter +1;
        else
            %increment the number
            pathOutfile = strrep(pathOutfile,['_',num2str(intCounter-1),'.mat'],['_',num2str(intCounter),'.mat']);
            intCounter = intCounter +1;
        end
    end
    if (intCounter ~= 2)
        fprintf('\nWarning: outfile exists, new file saved as _%i.mat\n',intCounter-1)
    end
    
    save(pathOutfile,'sig_matrix','actual_t_matrix','strSourceFile','lstPValues','lstExcludedElectrodes');
end

fprintf('\nDone\n\n')

clearvars -except 'matrixSubjectByElectrode_*' strFilename strPhase dirOutput

%% Analyze P-values across different interations
%
%This code cell is intended to investigate the relationships between the
%number of excluded subjects for a given electrode and the variability of
%the p-value for that electrode; the goal is to determine if the
%variability in excluded subjects have an impact on the variability we see
%in the permutated results

%user-defined vars
dirInfiles = '/nfs/erp-modaf/mc/MC-epoched-Feb11-Glenn/Cue/Headplots/Bootstrap/0-250ms x100 (for p-value investigation)/';

%analysis to do
boolStatsVSExcludedCount = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%begin algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get list of all files
lstInputFiles = dir([dirInfiles,'*.mat']);
lstInputFiles = struct2cell(lstInputFiles); %convert to struct
lstInputFiles = lstInputFiles(1,:);   %delete unwanted file information

%load in all data
matrixPValues = NaN(length(lstInputFiles),128);
matrixExcludedElectrodes = NaN(length(lstInputFiles),128);
for i = 1:length(lstInputFiles)
    load([dirInfiles,lstInputFiles{i}])
    matrixPValues(i,:) = lstPValues;
    matrixExcludedElectrodes(i,:) = lstExcludedElectrodes;
end

if boolStatsVSExcludedCount
    matrixStats = [];   %(excluded count, p-value stdev)
    for i = 0:max(matrixExcludedElectrodes(1,:))
        lstIndexes = find(ismember(matrixExcludedElectrodes(1,:),i));
        matrixData = matrixPValues(:,lstIndexes);
        matrixStats(end+1,1) = i;
        matrixStats(end,2) = std(matrixData(:));
        matrixStats(end,3) = var(matrixData(:));
    end
    idxNoData = find(isnan(matrixStats(:,2)));
    matrixStats(idxNoData,:) = [];
    plot(matrixStats(:,1),matrixStats(:,2)) %std
    figure
    plot(matrixStats(:,1),matrixStats(:,3)) %var
end

clear