% Temporal Order EEG: Pre-Processing to-eeg data
% Written by: Josh Phillips 1/2/2013
% Last Revision by: Evan Layher 1/11/2013


% MODIFY TRIGGERS


% This script uses Ben's strsearch function located in the path below
addpath('/nfs/to-eeg/code/functions/')

% This portion of the script modifies the triggers - specifically the
% trigger following the first crosshair of a trial. This trigger is modified
% to take into account performance. This is the trigger that
% will be used to define time point 0 for each epoch.
All_Events={EEG.event(:).type};
indx99=strsearch('99',All_Events);
temp = All_Events(indx99+7);
trigger_to_be_changed =indx99+1;
item_hit = strsearch('157',temp);
item_miss = strsearch('159',temp);
tempx = temp;
tempx(item_hit) = {'333'};
tempx(item_miss) = {'444'};
order_hit = strsearch('57',tempx);
order_miss = strsearch('59',tempx);

All_Events_Check = All_Events;
All_Events_Check(trigger_to_be_changed(item_hit)) = {'200'}; % Item hit triggers are relabeld 200
All_Events_Check(trigger_to_be_changed(item_miss)) = {'210'}; % Item miss triggers are relabeld 210
All_Events_Check(trigger_to_be_changed(order_hit)) = {'220'}; % Order hit triggers are relabeld 220
All_Events_Check(trigger_to_be_changed(order_miss)) = {'230'}; % Order miss triggers are relabeld 230

for i = 1:length(EEG.event)
EEG.event(i).type = cell2mat(All_Events_Check(i));
EEG.urevent(i).type = cell2mat(All_Events_Check(i));
end

disp('Well done modifying those triggers Evan and Josh! Now load in the channel locations.')

%% LOAD IN CHANNEL LOCATIONS
EEG.chanlocs=readlocs('/nfs/cmeteeg/data/neuroscan_64_cap_3_2_2011.ced'); % Load in Channel Location file
disp('The channel locations are loaded in. Now remove the EKG and EMG electrodes.')
%% REMOVE ELECTRODES THAT WERE NOT USED
EEG=pop_select(EEG, 'nochannel', {'EKG' 'EMG'}); % Remove electodes that were not used to record data
disp('EKG and EMG Electrodes are removed. Now downsample the data to 500hz.')
%% RESAMPLE AT 500hz
EEG=pop_resample(EEG, 500); % Downsample to 500hz
disp('Downsampling to 500hz completed. Now Re-reference data to the mastoids (Electrodes 33 and 43).')
%% REFERENCE DATA TO THE MASTOIDS
ref_chan = [33 43]; % Define reference channels (In our the case the mastoids)
EEG = pop_reref( EEG, ref_chan, 'keepref', 'on' ); % Re-reference to the ref_chan channels defined above (we are trying to take the average reference of the mastoids)
disp('Electrodes referenced to the mastoids (Electrodes 33 and 43). Now place a high pass filter of 0.5hz on the data.')
%% HIGH PASS FILTER AT 0.5HZ
EEG.data=eegfilt(EEG.data, EEG.srate, .5, []); % Use a 0.5hz HPF on the data
disp('The data has been high pass filtered at 0.5hz. Now Epoch the data')
%% EPOCH DATA INTO 200 TRIALS
EEG=pop_epoch(EEG, {'200' '210' '220' '230'}, [-1 17], 'epochinfo', 'yes'); % Epoch the data (should be 200 if the subject answered each trial)
EEG=pop_rmbase(EEG, [-1000 0]); % Remove the baseline (always use an input of millisceonds for the pop_rmbase function, while the function pop_epoch uses seconds...be aware)
disp('The Data should be epoched into 200 trials. Make sure that all 200 exist.If so, SAVE DATA. If not, TROUBLESHOOT. After saving, hand reject the data.')
%%
%HAND REJECT

% Save your hand rejected dataset before starting ICA.
%% ICA
EEG = pop_runica(EEG, 'icatype','runica','dataset',1,'options',{'extended' 1},'chanind',[1:64] );
disp('ICA is finished. SAVE YOUR DATA NOW!')
%% VIEW ICA Components - 
pop_prop(EEG,0,[1:15])

disp('The first 15 components are ready to be viewed. Take out the eye blink components and input those numbers into the next function.')
%% Remove Components and then save

EEG = pop_subcomp(EEG,[1 2 5 8],0); 
disp('ICA components removed. Save your data. Now do an automated rejection')
%% Automated Rejection

[EEG EEG.rejected_indexes] = pop_eegthresh(EEG,1,[1:64],-100,100,EEG.xmin,EEG.xmax,1,1);
disp('Do one last hand rejection then SAVE THE FINAL DRAFT! Congratulations! You finished your dataset.')