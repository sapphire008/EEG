%separate task epoches
%copy from EEG file
N=432; %number of events in the first task

EEG_WM=EEG;
EEG_4POP=EEG;

%event field
EEG_WM.event=EEG.event(1:N);
EEG_4POP.event=EEG.event((N+1):end);
for n = 1:length(EEG_4POP.event)
    EEG_4POP.event(n).urevent=EEG_4POP.event(n).urevent-N;
end

%data field
task_separation_timing=mean([EEG.urevent(N+1).latency,EEG.urevent(N).latency]);
EEG_WM.data=EEG.data(:,1:floor(task_separation_timing));
EEG_4POP.data=EEG.data(:,ceil(task_separation_timing):end);

%urevent field
EEG_WM.urevent=EEG.urevent(1:N);
EEG_4POP.urevent=EEG.urevent((N+1):end);

for n = 1:length(EEG_4POP.urevent)
    EEG_4POP.urevent(n).latency=EEG_4POP.urevent(n).latency-ceil(task_separation_timing)+1;
    EEG_4POP.event(n).latency=EEG_4POP.urevent(n).latency;
end

%pnts field
EEG_WM.pnts=size(EEG_WM.data,2);
EEG_4POP.pnts=size(EEG_4POP.data,2);

%xmax field
EEG_WM.xmax=EEG_WM.pnts/EEG_WM.srate;
EEG_4POP.xmax=EEG_4POP.pnts/EEG_4POP.srate;


EEG_WM = pop_saveset(EEG_WM,'filepath',...
    '/nfs/jong_exp/EEG_PFC/subjects/PFC101_031813/WM/',...
    'filename',[params.subjectID,'-pre-epoched.set']);
EEG_4POP = pop_saveset(EEG_4POP,'filepath',...
    '/nfs/jong_exp/EEG_PFC/subjects/PFC101_031813/4POP/',...
    'filename',[params.subjectID,'-pre-epoched.set']);

%% Remove WM block 3 and 4
M=216;%end fo block 2
EEG_WM=EEG;
task_separation_timing=mean([EEG.urevent(M+1).latency,EEG.urevent(M).latency]);


EEG_WM.event=EEG.event(1:M);
EEG_WM.urevent=EEG.urevent(1:M);
EEG_WM.data=EEG.data(:,1:floor(task_separation_timing));
EEG_WM.pnts=size(EEG_WM.data,2);
EEG_WM.xmax=EEG_WM.pnts/EEG_WM.srate;

EEG_WM = pop_saveset(EEG_WM,'filepath',...
    '/nfs/jong_exp/EEG_PFC/subjects/PFC101_031813/WM/',...
    'filename',['PFC101_031813','-pre-epoched.set']);

%% correct 4POP
N=432;
for n = 1:length(EEG.event)
    EEG.event(n).urevent=EEG.event(n).urevent-N;
end

EEG = pop_saveset(EEG,'filepath',...
    '/nfs/jong_exp/EEG_PFC/subjects/PFC101_031813/4POP/',...
    'filename',['PFC101_031813','-pre-epoched.set']);














