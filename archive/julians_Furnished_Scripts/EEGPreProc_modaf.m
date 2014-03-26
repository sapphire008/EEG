% set the paths and variables so they match the group

 eval(['load ',elecs_dir,'/bad_elecs_',group,'.mat'])
     
if s >= 100
    zeroes = ''
else

    if s >= 10
        zeroes = '0'

    else

        zeroes = '00'
         
    end

end

%if dose1
if isequal(dose,'dose-1') || isequal(dose,'dose1')
    bad_elecs = bad_elecs{s,2} 
end
%if dose2
if isequal(dose,'dose2')
    bad_elecs = bad_elecs{s,3} 
end

if isequal(dose,'dose3')
    bad_elecs = bad_elecs{s,4} 
end

bad_elecs = bad_elecs(bad_elecs <= 128)

%DS2
if isequal(DS,'DS2')

cd(preprocdir)
mkdir([num2str(group),num2str(zeroes),num2str(s)])
cd([num2str(group),num2str(zeroes),num2str(s)])


%load CNT file
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

%mc
if isequal(group,'mc')
EEG = pop_loadcnt([raw_eegdir,num2str(group),num2str(zeroes),num2str(s),'-pop-re-marked-400ms-baseline-',dose,'.cnt'], 'dataformat', 'int32', 'keystroke', 'on');
end

%ms
if isequal(group,'ms')
EEG = pop_loadcnt([raw_eegdir,'/ms/',num2str(group),num2str(zeroes),num2str(s),'-',dose,'-pop-re-marked-400ms-baseline.cnt'], 'dataformat', 'int32', 'keystroke', 'on');
end

%pc
if isequal(group,'pc')
EEG = pop_loadcnt([raw_eegdir,num2str(group),num2str(zeroes),num2str(s),'-pop-re-marked-400ms-baseline-',dose,'.cnt'], 'dataformat', 'int32', 'keystroke', 'on');
end

% load electrode locations
EEG.chanlocs=readlocs([elecs_dir,'/Oct-06-channel-locs.ced']);

%remove extra electrodes
    x= length(EEG.chanlocs)
    if x>128
    EEG = pop_select( EEG, 'nochannel', [129:x]); 
    end
    
%reject electrodes    
EEG = pop_select( EEG, 'nochannel', bad_elecs);

% downsample data to make the file smaller
EEG = pop_resample(EEG,250);

% high-pass filter at 0.5 Hz
EEG = pop_eegfilt( EEG, 0.5, 0, [], [0]);

%normal delay
if isequal(phase,'Cue')
EEG = pop_epoch( EEG, {'21'  '11'}, [-0.4 1.7], 'epochinfo', 'yes');
end

%probe-phase 
if isequal(phase,'Probe')

original = cell2mat({EEG.event.latency})/EEG.srate

for i = 1:2:length(EEG.event)

    old_lat = EEG.event(i).latency
    
    addition = 1.5*EEG.srate
    
    new_lat = old_lat+addition

    EEG.event(i).latency = new_lat
    
end

new_latencies = cell2mat({EEG.event.latency})/EEG.srate

check = [new_latencies - original]

EEG = pop_epoch( EEG, {'20' '21' '10' '11'}, [-0.4 1.2], 'epochinfo', 'yes');

end

% average reference (specific to Carter Lab Neuroscan caps)
EEG = pop_reref( EEG, [], 'refstate',0);

% baseline subtract
EEG = pop_rmbase( EEG, [-400 0]);

% save the epoched data set
EEG = pop_saveset( EEG,  'filename', [num2str(group),num2str(zeroes),num2str(s),'-',dose,'-epoched.set'], 'filepath', [preprocdir,num2str(group),num2str(zeroes),num2str(s)]);

clear EEG;

end


%DS1
if isequal(DS,'DS1')

cd(preprocdir)
mkdir([num2str(group),num2str(zeroes),num2str(s)])
cd([num2str(group),num2str(zeroes),num2str(s)])

%NOTES- since the DS1 subjects cannot be remarked, the raw .cnt files must be MANUALLY imported
%via EEGLAB. Dialog box will pop up.

EEG = pop_loadcnt()

% load electrode locations
EEG.chanlocs=readlocs([elecs_dir,'/Oct-06-channel-locs.ced']);

x= length(EEG.chanlocs)
if x>128
    EEG = pop_select( EEG, 'nochannel', [129:x]); 
end

EEG = pop_select( EEG, 'nochannel', bad_elecs);


EEG = pop_resample(EEG,250);

% % high-pass filter at 0.5 Hz
EEG = pop_eegfilt( EEG, 0.5, 0, [], [0]);
EEG = pop_saveset( EEG,  'filename', [num2str(group),num2str(zeroes),num2str(s),'-',dose,'-unepoched.set'], 'filepath', [preprocdir,num2str(group),num2str(zeroes),num2str(s)]);

clear EEG

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%% for ds1, remark using the ev2 file %%%%%%
%%%%%%%%%%%%%     GLENN VERSION      %%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

EEG = pop_loadset( 'filename', [num2str(group),num2str(zeroes),num2str(s),'-',dose,'-unepoched.set'], 'filepath', [preprocdir,num2str(group),num2str(zeroes),num2str(s)]);
% 
clearvars ev2 data textdata colheaders ev2dir

ev2dir = '/nfs/erp-modaf/evt-files/mc_evt/csvs/';
cd(ev2dir)

importdata([ev2dir,group,zeroes,num2str(s),'-',dose,'-ev2.csv'],',',1)

ev2 = data;    

N = length(ev2);
    new_lat = zeros(1,N);

    for i = 1:N
        new_lat(i) = round((ev2(i,12)*250)/1000); % convert to samples
        %new_lat(i) = round((ev2(i,7)*250)/1000); % convert to samples
    end

clearvars ev2 data textdata colheaders ev2dir    
    
ev2dir = '/nfs/erp-modaf/evt-files/evt-files/';
cd(ev2dir)

if isequal(dose,'dose1')
    dose_num = '1st'
end
if isequal(dose,'dose2')
    dose_num = '2nd'
end
importdata(eval([ev2dir,group,zeroes,num2str(s),'-',dose_num,'-pop-stimulus-corrected.ev2']),',',1)



ev2 = data;     
N = length(ev2);
    for i = 1:N
        EEG.event(1,i).latency = new_lat(i);
        %EEG.event(1,i).type = num2str(ev2(i,8)); % 2nd col is eventtype
        EEG.event(1,i).type = ev2(i,2); % 2nd col is eventtype
    end
 
    
    
 
EEG = pop_saveset( EEG,  'filename', [num2str(group),num2str(zeroes),num2str(s),'-',dose,'_fixed.set'], 'filepath', [preprocdir,num2str(group),num2str(zeroes),num2str(s)], 'savemode', 'onefile');
clear EEG 

end


% Artifact Rejection

cd(preprocdir)
%load the epoched data
EEG = pop_loadset( 'filename', [num2str(group),num2str(zeroes),num2str(s),'-',dose,'-epoched.set'], 'filepath', [preprocdir,num2str(group),num2str(zeroes),num2str(s)]);


% reject epochs based on probability
[EEG, locthresh, globthresh, nrej] = pop_jointprob(EEG,1,[1:EEG.nbchan] ,3,5,1,0,0);

%Reject marked epochs and save the dataset
cd([num2str(group),num2str(zeroes),num2str(s)])
EEG = eeg_rejsuperpose( EEG, 1, 1, 1, 1, 1, 1, 1, 1);
EEG = pop_rejepoch( EEG, find(EEG.reject.rejglobal), 0);

EEG = pop_eegthresh(EEG,1,[1:EEG.nbchan],-150,150,EEG.xmin,EEG.xmax,1,0);

EEG = pop_saveset(EEG, 'filename', [num2str(group),num2str(zeroes),num2str(s),'-',dose,'-PostR.set'], 'filepath', [preprocdir,num2str(group),num2str(zeroes),num2str(s)]);

clear EEG locthresh globthresh nrej bad_elecs