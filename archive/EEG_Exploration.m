
%base work directory, where the preprocessing is saved
work_dir='/nfs/jong_exp/EEG_4POP/subjects/pilot_02012013/';
load([work_dir 'workspace.mat']);
addpath(genpath('/home/cui/eeg_4pop/test_1/scripts'));
addpath(genpath('/home/cui/scripts/EEGLABv11'));
%% Exploration 1: Spectrogram
load([work_dir 'Data/ERSP_pilot_991_30-80Hz.mat']);
load([work_dir 'Data/ITC_pilot_991_30-80Hz.mat']);
load([params.Save_dir,'workspace.mat']);
cond_type = {'incong', 'cong', 'contrast'};%a list of condition types
PFC_type = Region(:,1);%a list of PFC regions
clear i j k;

% Phase x Condition x Brain Region x Freq x Time
for i = 1:numPhase%Phase level
    %load in necessary parameters for each Phase
    times_out=TFparams.(Phase{i,1}).times_out;
    freqs_out=TFparams.(Phase{i,1}).freqs_out;
    for j = 1:length(cond_type) %Condition level
        %find index available electrodes (ones having freq x time data)
        ERSP_avail_elecs=find(~cellfun(@isempty,ERSP.(Phase{i,1}).(cond_type{1,j}).lead));
        for k=1:length(PFC_type) %Brain Region level
            %get electrodes that are both available and 
            get_elecs=intersect(Dict.Region.(PFC_type{k}),ERSP_avail_elecs);
            %store the freq x time data to the spectogram
            ERSP_Spectrogram.(Phase{i,1}).(cond_type{j}).(PFC_type{k})=...
                ERSP.(Phase{i,1}).(cond_type{j}).lead(get_elecs);
            %average freq x time spectogram across all electrodes at a
            %region
            ERSP_Spectrogram_elecs_avg.(Phase{i,1}).(cond_type{j}).(...
                PFC_type{k})=EEG_meanSpect(ERSP_Spectrogram.(...
                Phase{i,1}).(cond_type{j}).(PFC_type{k}));
            %generating a spectogram image
            plot_name=[Phase{i,1},'-'....
                cond_type{j},'-',PFC_type{k}];%title of the plot
            h = tftopo(ERSP_Spectrogram_elecs_avg.(Phase{i,1}).(...
                cond_type{j}).(PFC_type{k}),times_out,freqs_out,...
                'title',plot_name,'verbose','off');
            caxis([-1.5 1.5]);
            set(gca,'YTick',[freqs_out(1):5:freqs_out(end)]);
            colorbar
            figure;
            saveas(1, [params.Save_dir,'Spectrogram/',...
                plot_name,'.tif'],'tiff');
            close figure;
        end 
    end
end

save([params.Save_dir,'Spectrogram/gamma/','ERSP_Spectrogram.mat'],...
    'ERSP_Spectrogram','ERSP_Spectrogram_elecs_avg');
%% Exploration 2: HeadPlot (on hold)
mkdir([params.Save_dir,'HeadPlot']);
load([params.Save_dir,'Data/RawDataSet']);
load([params.Save_dir,'Data/ERSP_Spectrogram.mat']);
% set scale
scale_limits = [-1.5 1.5]
%create headplot spline
clear i;
for i = 1:numPhase
    headplot('setup',RawDataSet.(Phase{i,1}).chanlocs,...
        [Phase{i,1},'headplot.spl']);
    
    figure; 
    headplot(EEG.data,'headplot.spl','view',[0 90],...
        'maplimits',scale_limits,'cbar',0,'electrodes','off');
    %saveas([params.Save_dir,'HeadPlot/',Phase{i,1},'-HeadPlot.tif'],...
        %'tiff');
end


EEG = pop_loadset('/nfs/erp-stroop/Glenn-Mar11/Cued_Stroop/sc-generic_set.set')

% remove electrode locations by matching electrode number with EEG.chanlocs
%.type

bad_elecs = find(isnan(cue_mean_gamma_struct(:,1))')

elecs = cellfun(@str2num,{EEG.chanlocs.type})

remove_elecs = []
for b=1:length(bad_elecs)
    if ~isempty(find(elecs == bad_elecs(b)))
        remove_elecs(end+1) = find(elecs == bad_elecs(b))
    end
end

EEG = pop_select( EEG, 'nochannel', remove_elecs);

cue_mean_gamma_struct(bad_elecs,:) = []

EEG.data = mean(cue_mean_gamma_struct,2)

headplot('setup', EEG.chanlocs, 'headplot.spl');

figure; headplot(EEG.data,'headplot.spl','view',[0 90],'maplimits',scale_limits,'cbar',0,'electrodes','off')    




