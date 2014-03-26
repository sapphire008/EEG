%function EEG_GroupAnalysis(ERSP)
%% PART I: Set Processing Parameters
% reorganize data into <subjects x elecs x freqeuncy x time> 4D matrix
% length elecs = 63 consistently

% Specify processing directory
base_dir='/nfs/jong_exp/EEG_PFC/';%source of data
save_dir=[base_dir,'GroupAnalysis/'];%save directory of processed data
dataType='ERSP';
% Get group level Data information
%list_dir=dir([base_dir,'PFC*']);
%subjects=arrayfun(@(x) x.name, list_dir,'UniformOutput',false);
%list subjects
subj_list={...
    'PFC101_031813','PFC102_032013',...
    'PFC103_040313','PFC104_041213',...
    'PFC200_030813','PFC201_032813',...
    'PFC202_032213','PFC203_040813',...
    'PFC204_041713'};
template_subj = 'PFC104_041213';%used to get data structure and information
subj_group = {'SZ',200:299;'C',100:199};
flag.SeparateGroup = 0;
%list tasks to be processed
Task_list={'4POP','WM'};

% ======The following does not necessarily need to be modified ======
%get data info
%note that if the field names contain numbers, it will be converted into
%English spelling of the number
DataInfo=EEG_GetDataInfo(base_dir,Task_list,template_subj);
DataInfo.subjects = subj_list;
DataInfo.subj_group=subj_group;
DataInfo.subj_vect=EEG_groupSubjects(subj_list,subj_group,'PFC');
%to use in structure, we have to convert number in the task name to English
%spellings of the number
tmp_task=cellfun(@Number2Word,Task_list,'UniformOutput',0);
%dimensions of the data matrix
dim_name='<subject x elecs x frequency x time>';
save([save_dir,'Group_analysis_workspace.mat']);

load('/nfs/jong_exp/EEG_PFC/scripts/modules/defaultParams.mat','Dict');
Region_list = fieldnames(Dict.Region64);

clear subj_list subj_group template_subj;

%% PART II: Create and Save Group Database
%for each task(stored in separate folders)
clear tk;
for tk = 1:length(Task_list)
    %make a folder for each task
    mkdir([save_dir,Task_list{tk}]);
    clear Elecs Phases Conditions;
    %unwinding DataInfo
    Elecs=DataInfo.(tmp_task{tk}).Elecs;
    Phases=DataInfo.(tmp_task{tk}).Phases;
    Conditions=DataInfo.(tmp_task{tk}).Conditions;
    %log skipped subjects
    skipped_subj.(tmp_task{tk}) = {};
    
    %for each subject group (stored in separate folders under task folder
    for g = 1:size(DataInfo.subj_group,1)
        %flag whether or not to use group separation
        if flag.SeparateGroup
            subjects = DataInfo.subjects(DataInfo.subj_vect.(...
                DataInfo.subj_group{g}));
            group_tag = DataInfo.subj_group{g};
        else
            subjects = DataInfo.subjects;
            group_tag = 'All';
        end
        %make folders for each group
        mkdir([save_dir,Task_list{tk},'/',group_tag,...
            '/Spectrogram_gamma']);
        
        %for each phase(stored in separate files)
        clear n;
        for n = 1:length(Phases)
            %for each conditions (stored in separated fields)
            clear m;
            Group_RgnAVG_Data = struct();
            Group_SubjAVG_Data = struct();
            for m = 1:length(Conditions)
                %place hold data structure for all conditions
                GroupDataSet.(Conditions{m})=zeros(...
                    length(subjects),...
                    Elecs,...
                    DataInfo.(tmp_task{tk}).(Phases{n}).Freqs,...
                    DataInfo.(tmp_task{tk}).(Phases{n}).Times);
                %load and concatenate each subject's data
                
                clear s;
                for s = 1:length(subjects)
                    get_file_name=dir([base_dir,'subjects/',...
                        subjects{s},'/',Task_list{tk},'/Data/',...
                        dataType,'_PFC*.mat']);
                    
                    %load data if exists
                    if isempty(get_file_name)
                        skipped_subj.(tmp_task{tk}){end+1}=subjects{s};
                        clear get_file_name;
                        continue;%skip
                    else%if exist proceed to load the data
                        tmp_data=load([base_dir,'subjects/',subjects{s},...
                            '/',Task_list{tk},'/Data/',...
                            get_file_name.name]);%e.g. tmp_data.ERSP
                        clear get_file_name;
                    end
                    
                    %store data to GroupDataSet structure
                    %fill empty cell with nan
                    tmp_nan_filled=FillEmptyCell(tmp_data.(dataType).(...
                        Phases{n}).(Conditions{m}).lead);
                    clear tmp_data;
                    %store data
                    GroupDataSet.(Conditions{m})(s,:,:,:)=shiftdim(cat(...
                        3,tmp_nan_filled{:}),2);
                    clear tmp_nan_filled;
                end%end subject

                %for each brain region, generate a plot
                for re = 1:length(Region_list)
                    %take mean over regions defined by electrodes
                    Group_RgnAVG_Data.(Conditions{m}).(...
                        Region_list{re}) = squeeze(nanmean(...
                        GroupDataSet.(Conditions{m})(:,Dict.Region64.(...
                        Region_list{re}),:,:),2));
                    %take mean over subjects
                    Group_SubjAVG_Data.(Conditions{m}).(...
                        Region_list{re}) = squeeze(nanmean(...
                        Group_RgnAVG_Data.(Conditions{m}).(...
                        Region_list{re}),1));
                    
                    % Plot
                    %title of the plot
                    plot_name=[Phases{n},'-'....
                        Conditions{m},'-',Region_list{re}];
                    % plot spectrogram
                    EEG_PlotSpectrogram(...
                        Group_SubjAVG_Data.(...
                        Conditions{m}).(Region_list{re}),...
                        DataInfo.(tmp_task{tk}).(Phases{n}).times_out,...
                        DataInfo.(tmp_task{tk}).(Phases{n}).freqs_out,...
                        [save_dir,Task_list{tk},'/',group_tag,...
                        '/Spectrogram_gamma/'],plot_name)
                    clear plot_name;
                end%end region
            end%end condition
            
            skipped_subj.(tmp_task{tk})=unique(skipped_subj.(tmp_task{tk}));
            %save data of current phase of current group of current task
            DataInfo.(tmp_task{tk}).(Phases{n}).(...
                group_tag).data_dir{1} = ...
                [save_dir,Task_list{tk},'/',group_tag,...
                '/EEG_PFC_Group_',dataType,'_',...
                Phases{n},'-RAW.mat'];
            save(DataInfo.(tmp_task{tk}).(Phases{n}).(...
                group_tag).data_dir{1},...
                'GroupDataSet','dim_name','subjects','skipped_subj',...
                '-v7.3');
            DataInfo.(tmp_task{tk}).(Phases{n}).(...
                group_tag).data_dir{2} = ...
                [save_dir,Task_list{tk},'/',group_tag,...
                '/EEG_PFC_Group_',dataType,'_',...
                Phases{n},'-RegionAVG.mat'];
            save(DataInfo.(tmp_task{tk}).(Phases{n}).(...
                group_tag).data_dir{2},...
                'Group_RgnAVG_Data');
            DataInfo.(tmp_task{tk}).(Phases{n}).(...
                group_tag).data_dir{3} = ...
                [save_dir,Task_list{tk},'/',group_tag,...
                '/EEG_PFC_Group_',dataType,'_',...
                Phases{n},'-SubjAVG.mat'];
            save(DataInfo.(tmp_task{tk}).(Phases{n}).(...
                group_tag).data_dir{3},...
                'Group_SubjAVG_Data');
            save([save_dir,'Group_analysis_workspace.mat'],...
                'DataInfo','skipped_subj','-append');
            clearvars Group_*;
        end%end phase
        if strcmpi(group_tag,'All')
            break;%termoinate loop if running all subjects
        end
    end%end group
end%end task

%% PART III: Permutation Test
workspace_dir = '/nfs/jong_exp/EEG_PFC/GroupAnalysis/ERSP_Group_analysis_workspace.mat';
%add path of GPU function
addpath /nfs/pkg64/gpu/
%load dictionary
load('/nfs/jong_exp/EEG_PFC/scripts/modules/defaultParams.mat','Dict');
%load data information
load(workspace_dir);
% Set up comparison conditionse
%use of table: 'Within' each condition, compare 'Between' conditions
%for instance, 'Within' incongruent conditions, compare SZ and C
%or, 'Within' SZ subject Group, compare incong and cong conditions
Comparison_table = {...
    'comparison_type', 'Within',                      'Between',           'ttest';...
    'Between_Groups',  {'incong','cong','contrast'},  {'C','SZ'},          'pooled';...
    'Within_Groups',   {'SZ','C'},                    {'incong','cong'},   'paired'};
Comparison_struct = Table2Struct(Comparison_table,2);
Comparison_Fields = EEG_makeComparisonFields(Comparison_struct,...
    {{'SZ','C'},{'incong','cong','contrast'}});
Region_list = fieldnames(Dict.Region64);
dblAlpha = 0.05;% significance level

% Begin Permutation Test (DO NOT MODIFY THE FOLLOWING)
%for each task
for tk = 2:length(Task_list)
    mkdir([save_dir,Task_list{tk},'/Permutation_test']);
    Phases = DataInfo.(tmp_task{tk}).Phases;
    %for each phase
    for n = 1:length(Phases)
        comparison_types = fieldnames(Comparison_struct);
        %for each comparison type
        for x = 1:length(comparison_types)
            PermutationData = struct();
            mkdir([save_dir,Task_list{tk},'/Permutation_test/',...
                comparison_types{x},'/Spectrogram_gamma']);
            %load data for each group
            %Data (.Group  .Condition  .Region . <subjects x freqs x time>)
            Data.SZ = getfield(load(DataInfo.(tmp_task{tk}).(...
                Phases{n}).SZ.data_dir{2}),'Group_RgnAVG_Data');
            Data.C = getfield(load(DataInfo.(tmp_task{tk}).(...
                Phases{n}).C.data_dir{2}),'Group_RgnAVG_Data');
            
            Within_comparisons = ...
                (Comparison_struct.(comparison_types{x}).Within);
            
            tmp = Comparison_Fields.(comparison_types{x});
            %for each comparison
            for y = 1:size(tmp)
                data_struct_1 = Data.(tmp{y,1}{1}).(tmp{y,1}{2});
                data_struct_2 = Data.(tmp{y,2}{1}).(tmp{y,2}{2});
                %only run regions that both data struct has data on
                common_regions = intersect(fieldnames(data_struct_1),...
                    fieldnames(data_struct_2));
                %for each region
                for re = 1:length(common_regions)
                    data_mat_1 = data_struct_1.(common_regions{re});
                    data_mat_2 = data_struct_2.(common_regions{re});
                    cycles = min(...
                        [2^size(data_mat_1,1),2^size(data_mat_2,1),5E6]);
                    p_value = zeros(size(data_mat_1,2),size(data_mat_1,3));
                    t_stats = p_value;
                    %run permutation test on current region of current
                    %condition of comparison of current comparison type of
                    %curernt phase of current task, on each freq x time
                    for i = 1:size(data_mat_1,2)
                        for j = 1:size(data_mat_1,3)
                            [p_value(i,j),t_stats(i,j)] = ...
                                permutation_test(...
                                data_mat_1(:,i,j)',data_mat_2(:,i,j)',...
                                Comparison_struct.(...
                                comparison_types{x}).ttest,cycles);
                        end
                    end
                    %store p-value and t-stat map
                    PermutationData.cycles = cycles;
                    PermutationData.ttest = Comparison_struct.(...
                        comparison_types{x}).ttest;
                    PermutationData.p_value.(comparison_types{x}).(...
                        [tmp{y,1}{1},'_',tmp{y,1}{2},'_vs_',...
                        tmp{y,2}{1},'_',tmp{y,2}{2}]).(...
                        common_regions{re}) = p_value;
                    p_value_alpha = p_value;
                    p_value_alpha(p_value>dblAlpha)=1;
                    PermutationData.p_value_alpha.(comparison_types{x}).(...
                        [tmp{y,1}{1},'_',tmp{y,1}{2},'_vs_',...
                        tmp{y,2}{1},'_',tmp{y,2}{2}]).(...
                        common_regions{re}) = p_value_alpha;
                    
                    PermutationData.t_stats.(comparison_types{x}).(...
                        [tmp{y,1}{1},'_',tmp{y,1}{2},'_vs_',...
                        tmp{y,2}{1},'_',tmp{y,2}{2}]).(...
                        common_regions{re}) = t_stats;
                    t_stats_alpha = t_stats;
                    t_stats_alpha(p_value>dblAlpha) = 0;
                    
                    PermutationData.t_stats_alpha.(comparison_types{x}).(...
                        [tmp{y,1}{1},'_',tmp{y,1}{2},'_vs_',...
                        tmp{y,2}{1},'_',tmp{y,2}{2}]).(...
                        common_regions{re}) = t_stats_alpha;
                    
                    
                    %plot and save figures
                    %tmp_cell: {data_to_graph,'save_name_suffix',[color
                    %range], reverse_color_switch}
                    half_color = jet;
                    half_color = half_color(round(size(half_color,1)/2):end,:);
                    tmp_cell = {p_value,'-P',[0,1],1,half_color;...
                            p_value_alpha,'-P_ALPHA',[0,1],1,half_color;...
                            t_stats,'-T',[-5,5],0,[];...
                            t_stats_alpha,'-T_ALPHA',[-5,5],0,[]};
                    for kk =1:size(tmp_cell,1)
                        EEG_PlotSpectrogram(tmp_cell{kk,1},...
                            DataInfo.(tmp_task{tk}).(Phases{n}).times_out,...
                            DataInfo.(tmp_task{tk}).(Phases{n}).freqs_out,...
                            [save_dir,Task_list{tk},'/Permutation_test/',...
                            comparison_types{x},'/Spectrogram_gamma/'],...
                            [Phases{n},'-',[tmp{y,1}{1},'_',...
                            tmp{y,1}{2},'_vs_',tmp{y,2}{1},'_',tmp{y,2}{2}],...
                            '-',common_regions{re},tmp_cell{kk,2}],...
                            'ColorRange',tmp_cell{kk,3},...
                            'ReverseColor',tmp_cell{kk,4},...
                            'ColorMap',tmp_cell{kk,5});
                    end
                    
                    clearvars t_st* p_va* tmp_cell kk;
                end%end regions
            end%end comparisons
            
            save([save_dir,Task_list{tk},'/Permutation_test/',...
                comparison_types{x},'/',dataType,'-PermutationData',...
                '-',Phases{n},'.mat'],'PermutationData');
            clear PermutationData;
        end %end comparison types
    end%end phases
end%end task


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