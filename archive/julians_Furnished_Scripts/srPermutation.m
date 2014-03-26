%% bootstrap ERSP (adapted from Josh's version)
%This file is a modified version of CMET_bootstrap_corrected_between_gpu.m
%from the SC study, which was originally given to Julian Cheng by Josh
%Phillips

%LOAD INPUT 4D ERSP/ITC MAT FILES BEFORE PROCEEDING
%MAKE SURE THIS IS RUN ON A MACHINE WITH GPU HARDWARE

%Changelog:
%   2/5/2013:   changed behavior to handle complex values for ITC. added
%               functionality to automatically parse the input filename and
%               determine the data type and frequency range
%   2/6/2013:   reverted all changes done on 2/5/2013 that involves ABS()
%               of ITC values; this was incorrect

%% Perform permutations

addpath /home/phillips/Scripts/
addpath /nfs/pkg64/gpu

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

%user-defined vars

%select comparison group to run
%0 = Drug Red vs Green
%1 = Plac Red vs Green
%2 = Q money, i.e. Drug contrast vs Plac contrast
%3 = Other, makes use of lstCustomGroup
flagGroup = 2;
lstCustomGroup = [{'Drug_incong'} {'Drug_cong'}];
strGroup = 'Drug_incongVcong';  %used in filename

%select regions to run
%region_average_names = [{'left_PFC'}];
region_average_names = [ {'left_PFC'} {'mid_PFC'} {'right_PFC'} {'left_Central'} {'mid_Central'} {'right_Central'} {'left_Post'} {'mid_Post'} {'right_Post'}];
%region_average_names = [ {'all_PFC'} {'left_PFC'} {'mid_PFC'} {'right_PFC'} ];
%region_average_names = [{'mid_PFC'}];

%set permutation parameters
perm_number = 50000;
alpha = .05;

%select data type *Note: affects file names of output
%Update: if strFilename is present in workspace, will parse that instead.
%        If this is undesirable, set boolAutoParse_DataType = false;
%0 = ERSP
%1 = ITC
flagDataType = 1;
boolAutoParse_DataType = true;

%select frequency range *Note: affects file names of output
%Update: if strFilename is present in workspace, will parse that instead.
%        If this is undesirable, set boolAutoParse_Frequency = false;
%0 = 4-30Hz
%1 = 30-80Hz
flagFrequency = 0;
boolAutoParse_Frequency = true;

%set output directory
%Update: if strPhase is present in workspace, will use that to change between phases.
%        If this is undesirable, set boolAutoParse_DirOutput = false;
dirOutput = '/nfs/erp-modaf/mc/MC-epoched-Feb11-Glenn/Cue/New bootstrapping/Output/';
boolAutoParse_DirOutput = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%begin algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%auto-parse inputs from workspace variables
if boolAutoParse_DataType && exist('strFilename','var')
    strSplit = regexp(strFilename,'_','split');
    strDataType = strSplit{1};
    if strcmp(strDataType,'ERSP')
        flagDataType = 0;
    elseif strcmp(strDataType,'ITC')
        flagDataType = 1;
    else
        fprintf('Warning: failed to parse strFilename, using defined value for flagDataType\n')
    end
elseif boolAutoParse_DataType
    fprintf('Warning: strFilename is not present, using defined value for flagDataType\n')
end
if boolAutoParse_Frequency && exist('strFilename','var')
    strSplit = regexp(strFilename,'_','split');
    strTemp = strSplit{4};
    idxHzLabel = strfind(strTemp,'Hz');
    strFrequency = strTemp(1:idxHzLabel-1);
    if ~isempty(strfind(strFrequency,'4')) && ~isempty(strfind(strFrequency,'30'))
        flagFrequency = 0;
    elseif ~isempty(strfind(strFrequency,'30')) && ~isempty(strfind(strFrequency,'80'))
        flagFrequency = 1;
    else
        fprintf('Warning: failed to parse strFilename, using defined value for flagFrequency\n')
    end
elseif boolAutoParse_Frequency 
    fprintf('Warning: strFilename is not present, using defined value for flagFrequency\n')
end
if boolAutoParse_DirOutput && exist('strPhase','var')
    dirOutput = ['/nfs/erp-modaf/mc/MC-epoched-Feb11-Glenn/',strPhase,'/New bootstrapping/Output/'];
elseif boolAutoParse_DirOutput
    fprintf('Warning: strPhase is not present, using defined value for dirOutput\n')
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
    otherwise
        group = lstCustomGroup;
end

switch flagDataType
    case 0
        strDataType = 'ERSP';
    case 1
        strDataType = 'ITC';
    otherwise
        error('Unknown data type selection')
end

switch flagFrequency
    case 0
        strFrequency = '4-30Hz';
    case 1
        strFrequency = '30-80Hz';
    otherwise
        error('Unknown frequency selection')
end

for l = 1:length(region_average_names)
    
    eval(['region_average_elecs =', region_average_names{l},';'])

elecs = 1;

    % squeeze and nanmean across region to create a single electrode
    % (average)
    eval(['first_group_contrast_region = squeeze(nanmean(',group{1},'(:,region_average_elecs,:,:),2));'])
    eval(['second_group_contrast_region = squeeze(nanmean(',group{2},'(:,region_average_elecs,:,:),2));'])

    max_f = size(eval([group{1}]),3);
    max_t = size(eval([group{1}]),4);
    
    actual_t_matrix = zeros(max_f, max_t, elecs);
    sig_matrix = zeros(max_f, max_t, elecs);

    tic
for elec = 1:elecs
    
    for f = 1:max_f
      
        for t = 1:max_t
                      
            cond1_vect = first_group_contrast_region(:,f,t);
            cond2_vect = second_group_contrast_region(:,f,t);         
           
            combined_vector = [cond1_vect ; cond2_vect];

            fraction = bootstrap(combined_vector, length(cond1_vect), length(cond2_vect), perm_number);

            clear combined_vector group_* sequence

            sample_mean_difference = (mean(cond1_vect)-mean(cond2_vect));

            actual_t_matrix(f,t,elec) = sample_mean_difference;

            if fraction < alpha/2 || fraction > (1 - alpha/2)
                    sig_matrix(f,t,elec) = sample_mean_difference;
                else
                    sig_matrix(f,t,elec)= 0;
            end
                
            clear h p ci stats upper_t lower_t pseudo* cond1_vect cond2_vect combined_vector sorted*
       end
       disp(['Region: ', region_average_names{l},' | Freq: ',num2str(f),'/',num2str(size(eval([group{1}]),3)),' | Time: ',num2str(t),'/',num2str(size(eval([group{1}]),4))])
    end
                
end
    toc

    
if exist('strPhase','var')
    save([dirOutput,strPhase,'_',strDataType,'_',strFrequency,'_',strrep(region_average_names{l},'_','-'),'_',strGroup,'_bootstrap'], 'sig_matrix', 'actual_t_matrix','new_freqs','new_times')
else
    save([dirOutput,strDataType,'_',strFrequency,'_',strrep(region_average_names{l},'_','-'),'_',strGroup,'_bootstrap'], 'sig_matrix', 'actual_t_matrix','new_freqs','new_times')
end

end

fprintf('\nDone\n\n')

clearvars -except '*cong' '*incong' '*contrast' 'new*' strPhase strFilename

%% Plot
%close all;
% it appears that sig_matrix hold the bootstrapped data, and
% actual_t_matrix hold the un-bootstrapped data
elec = 1; %set =1 for average region
strTitle = 'myTitle';
%arrCScale = [-0.33 0.33];    %used to scale coloration of graph
%arrCScale = [-0.033 0.033];       %ITC only
arrCScale = [0 0.08];

freqs = new_freqs;
time_interval = new_times;

%select whether to view all t-stats (actual_t_matrix) or bootstrapped
%results (sig_matrix)
flagSwitch = 0; %0 = sig matrix, 1 = actual t matrix
if flagSwitch
    sjset = squeeze(actual_t_matrix(:,:,elec));
else
    sjset = squeeze(sig_matrix(:,:,elec));
end

figure
h = tftopo(sjset,time_interval,freqs,'title',strTitle,'verbose','off');
    caxis(arrCScale);
    colorbar;
%% Multiplot
%Sequentially plots data and saves as tiff
%Note: ITC plots only work if filename has "ITC" in it

close all
clear all

%vars
flagGraphs2Process = 0; %0 = both, 1 = sigma matrix only, 2 = actual t matrix only
intAverageElectrodes = 1;  %set to 1 for average region
dirRoot = '/nfs/erp-modaf/mc/MC-epoched-Feb11-Glenn/';
dirCue = 'Cue/New bootstrapping/';
dirProbe = 'Probe/New bootstrapping/';
lstToDo = [ ...  %this will process all files found within these folders for both cue and probe; if something doesn't exist it will warn and skip  
    %{'ERSP 4-30Hz/'} ...
    %{'ERSP 30-80Hz/'} ...
    {'ITC 4-30Hz/'} ...
%     {'ERSP 4-30Hz MC-only/'} ...
%     {'ITC 4-30Hz MC-only/'} ...
    ];
arrCScale = [-0.33 0.33];   %lower dB for 128 channel cap
arrCScaleITC = [0 0.08];    %ITC only, old; from Glenn's script, but not used anymore
%arrCScaleITC = [-0.033 0.033]; %ITC only, new; for results produced with new script

figure

for dirPhase = [{dirCue} {dirProbe}]
    for i = 1:length(lstToDo)
        %build path to target folder and search contents
        dirThisFolder = [dirRoot,cell2mat(dirPhase),lstToDo{i}];
        lstBootstrapResults = dir([dirThisFolder,'*.mat']);
        lstBootstrapResults = struct2cell(lstBootstrapResults); %convert to struct
        lstBootstrapResults = lstBootstrapResults(1,:);   %delete unwanted file information

        if isempty(lstBootstrapResults)
            %does not exist, skip
            warning('%s does not exist for dir %s',lstToDo{i},cell2mat(dirPhase))
            continue
        end

        if ~exist([dirThisFolder,'Figures'],'dir')
            %dir for figures is absent, so make it
            mkdir(dirThisFolder,'Figures')
        end
        if ~exist([dirThisFolder,'Figures/','Sig matrix'],'dir')
            %dir for sigma matrix is absent, so make it
            mkdir([dirThisFolder,'Figures'],'Sig matrix')
        end
        if ~exist([dirThisFolder,'Figures/','Actual t matrix'],'dir')
            %dir for actual t matrix is absent, so make it
            mkdir([dirThisFolder,'Figures'],'Actual t matrix')
        end
        
        %loop through all files
        for j = 1:length(lstBootstrapResults)
            load([dirThisFolder,lstBootstrapResults{j}]);
            
            strTitle = lstBootstrapResults{j}(1:end-14);
            matrixSJSetSigma = squeeze(sig_matrix(:,:,intAverageElectrodes));
            matrixSJSetActualT = squeeze(actual_t_matrix(:,:,intAverageElectrodes));
            
            %fix scale for ITC
            if ~isempty(strfind(strTitle,'ITC'))
                myCScale = arrCScaleITC;
            else
                myCScale = arrCScale;
            end
            
            if (flagGraphs2Process == 1) || (flagGraphs2Process == 0)
                h = tftopo(matrixSJSetSigma,new_times,new_freqs,'title',strrep(strTitle,'_',' '),'verbose','off');
                caxis(myCScale), colorbar;
                saveas(gcf,[dirThisFolder,'Figures/Sig matrix/',strTitle,'.tif'],'tiff')
                clf
            end
            if (flagGraphs2Process == 2) || (flagGraphs2Process == 0)
                h = tftopo(matrixSJSetActualT,new_times,new_freqs,'title',strrep(strTitle,'_',' '),'verbose','off');
                caxis(myCScale), colorbar;
                saveas(gcf,[dirThisFolder,'Figures/Actual t matrix/',strTitle,'.tif'],'tiff')
                clf
            end
            
            clear actual_t_matrix sig_matrix new_freqs new_times
        end
    end
end

close,clear all