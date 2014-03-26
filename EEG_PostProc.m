function EEG_PostProc(work_dir)
%EEG time frequency decomposition and wavelet analysis
%% Set Necessary Parameters for Post-Processing
%assume pwd as working dir if not specified
if nargin<1
    work_dir=pwd;
elseif strcmpi(work_dir(end),'/')
    work_dir=work_dir(1:end-1);
end
load([work_dir,'/workspace.mat']);
addpath(genpath(params.Script_dir));
EEG_addeegscripts;

%% Extract Conditions (RawDataSet.mat)
clear n;
for n = 1:numPhase
%load .set file from pre-processing
EEG = pop_loadset('filepath',[params.Save_dir Phase{n,1}, ...
    '/'],'filename',[params.subjectID,'-Analysis-Ready.set']);
%passing necessary EEG paramters to a structure of RawDataSet
%RawDataSet stores only necessary parameters from EEG for analysis
RawDataSet.(Phase{n,1})=CopyField(EEG,{'pnts','xmax','xmin','chaninfo',...
    'chanlocs','srate'});
%preparing folders for data analysis
mkdir([params.Save_dir Phase{n,1} '/Conditions/']);

%make vectors of each condition
[Vector_struct,numCond,cond_code,cond_name]=...
    EEG_makeVector(EEG,Task,Condition,Phase(n,:),...
    'acc',Accuracy,'rmempty',1); %#ok<USENS,ASGLU>
RawDataSet.(Phase{n,1}).idx=Vector_struct.idx;
clear Vector_struct;

% Separating into conditions for each phase
for cond=1:numCond
    clear EEG;
    %load preproc .set file
    EEG = pop_loadset('filepath',[params.Save_dir Phase{n,1},'/'],...
        'filename',[params.subjectID,'-Analysis-Ready.set']);
    %select trials with current condition
    EEG = pop_select( EEG, 'trial', ...
        RawDataSet.(Phase{n,1}).idx.(Condition{cond,1})); 
    EEG = pop_saveset(EEG,...
        'filepath',[params.Save_dir Phase{n,1} '/Conditions/'],...
        'filename',[params.subjectID,'-',Phase{n,1},'_',...
        Condition{cond,1},'_all.set']);
    %passing necessary EEG data to a structure
    RawDataSet.(Phase{n,1}).data.(Condition{cond,1})=EEG.data;
end
end
save([params.Save_dir,'workspace.mat'],...
    'numCond','cond_code','cond_name','-append');
%save a copy of the extracted conditions in a .mat format
mkdir([params.Save_dir,'Data/']);
save([params.Save_dir,'Data/RawDataSet.mat'],'RawDataSet');
%% Time-Frequency Analysis (ERSP-*.mat,ITC-*.mat)
%NEED TO MODIFY: use time intervals for modaf(?)
%                need to increase time resolution(?)
%                need to modify to loop phase (check)

%set parameters
if ~exist('RawDataSet','var')%in case RawDataSet is not present
    load([params.Save_dir,'Data/RawDataSet.mat']);
end
if ~exist('Dict','var')%in case Dict is not present
    load([params.Module_dir,'defaultParams.mat']);
end
freqs_in=Dict.Frequency.(params.bandInterest).frequency;
cycles=Dict.Frequency.(params.bandInterest).cycles;
clear n Dict taskSet;

for n=1:numPhase
    disp(['Current Phase: ',Phase{n,1}]);
    %Time Frequency Decomposition
    [ERSP.(Phase{n,1}),ITC.(Phase{n,1}),powbase,times_out,freqs_out]=...
        EEG_TimeFreq(RawDataSet.(Phase{n,1}),Phase(n,:),Condition,...
        params.subjectID,freqs_in,cycles,...
        'TimeStep',5,'contrast','on');
     
    %save the parameters
    TFparams.(Phase{n,1}).powbase=powbase;
    TFparams.(Phase{n,1}).times_out=times_out;
    TFparams.(Phase{n,1}).freqs_out=freqs_out;
end
save([params.Save_dir,'Data/ERSP_', params.subjectID,'_',...
    num2str(min(freqs_out)),'-',num2str(max(freqs_out)),'Hz','.mat'],...
    'ERSP');
save([params.Save_dir,'Data/ITC_', params.subjectID,'_',...
    num2str(min(freqs_out)),'-',num2str(max(freqs_out)),'Hz','.mat'],...
    'ITC');
save([params.Save_dir,'workspace.mat'],'TFparams','-append');

%% Exploration 1: Spectrogram
% load([params.Save_dir 'Data/ERSP_',params.subjectID,'_',...
%     num2str(min(freqs_out)),'-',num2str(max(freqs_out)),'Hz.mat']);
% load([params.Save_dir 'Data/ITC_',params.subjectID,'_',...
%     num2str(min(freqs_out)),'-',num2str(max(freqs_out)),'Hz.mat']);
load([params.Module_dir,'defaultParams.mat'],'Dict');
cond_type =Dict.Conditions;%a list of condition types
Region = fieldnames(Dict.(['Region',num2str(elecs.numChan)]));
clear c n k;
if ~exist([params.Save_dir,'Spectrogram_',params.bandInterest],'dir')
    mkdir([params.Save_dir,'Spectrogram_',params.bandInterest]);
end

% Phase x Condition x Brain Region x Freq x Time
for n = 1:numPhase%Phase level
    %load in necessary parameters for each Phase
    times_out=TFparams.(Phase{n,1}).times_out;
    freqs_out=TFparams.(Phase{n,1}).freqs_out;
    for c = 1:length(cond_type) %Condition level
        %find index available electrodes (ones having freq x time data)
        ERSP_avail_elecs=find(...
            ~cellfun(@isempty,ERSP.(Phase{n,1}).(cond_type{c}).lead));
        for k=1:length(Region) %Brain Region level
            %get electrodes that are both available and 
            get_elecs=intersect(Dict.(['Region',num2str(...
                elecs.numChan)]).(Region{k}),ERSP_avail_elecs);
            %store the freq x time data to the spectogram
            ERSP_Spectrogram.(Phase{n,1}).(cond_type{c}).(Region{k})=...
                ERSP.(Phase{n,1}).(cond_type{c}).lead(get_elecs);
            %average freq x time spectogram across all electrodes at a
            %region
            ERSP_Spectrogram_elecs_avg.(Phase{n,1}).(cond_type{c}).(...
                Region{k})=EEG_meanSpect(ERSP_Spectrogram.(...
                Phase{n,1}).(cond_type{c}).(Region{k}));
            %generating a spectogram image 
            plot_name=[Phase{n,1},'-'....
                 cond_type{c},'-',Region{k}];%title of the plot
            
            EEG_PlotSpectrogram(...
                ERSP_Spectrogram_elecs_avg.(Phase{n,1}).(...
                cond_type{c}).(Region{k}),times_out,freqs_out,...
                [params.Save_dir,'Spectrogram_',...
                 params.bandInterest,'/',plot_name,'.tif'],...
                'title',plot_name);
           
            %Save the selected electrode at current PhasexCondxRegion
            elecs.spect_elecs.(Phase{n,1}).(cond_type{c}).(Region{k})=...
                get_elecs;
        end 
    end
end

save([params.Save_dir,'Data/','ERSP_Spectrogram.mat'],...
    'ERSP_Spectrogram','ERSP_Spectrogram_elecs_avg');
save([params.Save_dir,'workspace.mat'], 'elecs','-append');
%% Exploration 2: HeadPlot (on hold)
% mkdir([params.Save_dir,'HeadPlot']);
% load([params.Save_dir,'Data/RawDataSet']);
% load([params.Save_dir,'Data/ERSP_Spectrogram.mat']);
% % set scale
% scale_limits = [-1.5 1.5];
% %create headplot spline
% clear n;
% for n = 1:numPhase
%     headplot('setup',RawDataSet.(Phase{n,1}).chanlocs,...
%         [Phase{n,1},'headplot.spl']);
%     
%     figure; 
%     headplot(EEG.data,'headplot.spl','view',[0 90],...
%         'maplimits',scale_limits,'cbar',0,'electrodes','off');
%     %saveas([params.Save_dir,'HeadPlot/',Phase{n,1},'-HeadPlot.tif'],...
%         %'tiff');
% end
% 
% 
% EEG = pop_loadset('/nfs/erp-stroop/Glenn-Mar11/Cued_Stroop/sc-generic_set.set');
% 
% % remove electrode locations by matching electrode number with EEG.chanlocs
% %.type
% 
% bad_elecs = find(isnan(cue_mean_gamma_struct(:,1))')
% 
% elecs = cellfun(@str2num,{EEG.chanlocs.type})
% 
% remove_elecs = []
% for b=1:length(bad_elecs)
%     if ~isempty(find(elecs == bad_elecs(b)))
%         remove_elecs(end+1) = find(elecs == bad_elecs(b))
%     end
% end
% 
% EEG = pop_select( EEG, 'nochannel', remove_elecs);
% 
% cue_mean_gamma_struct(bad_elecs,:) = []
% 
% EEG.data = mean(cue_mean_gamma_struct,2)
% 
% headplot('setup', EEG.chanlocs, 'headplot.spl');
% 
% figure; headplot(EEG.data,'headplot.spl','view',[0 90],'maplimits',scale_limits,'cbar',0,'electrodes','off')    




