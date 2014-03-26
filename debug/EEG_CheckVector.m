function EEG_CheckVector(EEG_dir,RawDataSet_dir,workspace_dir)
EEG_dir='/nfs/jong_exp/EEG_PFC/subjects/PFC103_040313/WM/Probe/PFC103_040313-Analysis-Ready.set';
RawDataSet_dir='/nfs/jong_exp/EEG_PFC/subjects/PFC103_040313/WM/Data/RawDataSet.mat';
workspace_dir='/nfs/jong_exp/EEG_PFC/subjects/PFC103_040313/WM/workspace.mat';

current_phase='Probe';
switch current_phase
    case {'Cue'}
        lookup_dir=1;
    case {'Probe'}
        lookup_dir=-2;
end

EEG=pop_loadset(EEG_dir);
load(RawDataSet_dir);

%get event time series first, grouped by trial/epoch
clear epoch_eventtype;
epoch_eventtype=zeros(length(EEG.epoch),5);
for ep = 1:length(EEG.epoch)
    epoch_eventtype(ep,[1:2]+0*length(EEG.epoch(ep).eventtype))=cell2mat(EEG.epoch(ep).eventtype);
    epoch_eventtype(ep,[2:3]+1*length(EEG.epoch(ep).eventtype))=cell2mat(EEG.epoch(ep).eventurevent);
    epoch_eventtype(ep,1+length(EEG.epoch(ep).eventtype))=EEG.urevent(epoch_eventtype(ep,5)+lookup_dir).type;
end
if strcmpi(current_phase,'Probe')
    tmp=epoch_eventtype;
    clear epoch_eventtype;
    epoch_eventtype(:,1)=tmp(:,3);
    epoch_eventtype(:,2:3)=tmp(:,1:2);
end
%load(workspace_dir);
%quick working memory
clear ep Vector;
Vector=struct('OneFace',[],'TwoFace',[]);
for ep = 1:size(epoch_eventtype,1)
    if epoch_eventtype(ep,3)==14
        switch epoch_eventtype(ep,1)
            case {31,32}%oneface
                Vector.OneFace(end+1)=ep;
            case {33 34}%twoface
                Vector.TwoFace(end+1)=ep;
        end
    end
end

%compare made vector from raw data to the vectors just made
try
    compare_logical(1)=(sum(sort(Vector.OneFace,1,'ascend')-sort(RawDataSet.(current_phase).idx.OneFace,1,'ascend'))==0);
    compare_logical(2)=(sum(sort(Vector.TwoFace,1,'ascend')-sort(RawDataSet.(current_phase).idx.TwoFace,1,'ascend'))==0);
catch
    error('Vector is not made correctly');
end
if ~compare_logical(1)
    error('Vector 1 is not made correctly');
elseif ~compare_logical(2)
    error('Vector 2 is not made correctly');
end


end


function flag_issubset=Check_issubset(set_A,set_B)
%check if set_A is the subset of set_B
flag_issubset=isempty(intersect(set_A,setxor(set_A,set_B)));
end

  