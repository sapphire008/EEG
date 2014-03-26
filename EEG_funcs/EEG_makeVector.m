function [Vector_struct,numCond,cond_code,cond_name]=EEG_makeVector(...
    EEG,Task,Condition,phase_in,varargin)
%calculate target events for conditions that we want to separate into
%4POP: {6 5}
%WM: {[31 32],[33 34],[35 36]}
% Optional paramters:
%   'acc': provided Accuracy information, the vector will only includes the
%          trials with correct response
%   'rmempty': ['0'|'1'], remove empty fields

% flag: if needs to exclude wrong
IND_acc=find(strcmpi({'acc'},varargin),1);
if ~isempty(IND_acc)
    Accuracy=varargin{IND_acc+1};
    ACC_target=FindACC(Accuracy);
end
%flag: if needs to remove empty field
IND_rm=find(strcmpi({'rmempty'},varargin),1);
if ~isempty(IND_rm)
    rm_empty=varargin{IND_rm+1};
else
    rm_empty=0;%default: do not remove empty field
end
    
%find conditions of interest
[cond_code,cond_name,numCond]=FindCondEvent(Task,Condition,phase_in);

%place holding
for nk=1:numCond
    Vector_struct.idx.(Condition{nk,1})=[];
end

% examining each epoch and extract condition and accuracy info
for ep = 1:1:length(EEG.epoch)
     %current_condition is the event code formatted as double
     cond_ind=FindCurrentCondition(...
         EEG.epoch(ep),EEG.urevent,cond_code,1);%cue condition
    if ~isempty(IND_acc)
        acc_ind=FindCurrentCondition(...
            EEG.epoch(ep),EEG.urevent,Accuracy(:,2),-1);%accuracy
        if Accuracy{acc_ind,2}~=ACC_target
            continue;%skip this epoch if current ACC does not match
        end
    end
    %storing what epoch #s belong to this condition
    Vector_struct.idx.(cond_name{cond_ind})(end+1)=ep;
end
if rm_empty%remove empty fields
    [Vector_struct.idx,emptyIND]=RemoveEmptyField(Vector_struct.idx);
    %recalculate number of conditions
    numCond=length(fieldnames(Vector_struct.idx));
end
end


%%
function [cond,cond_name,numCond]=FindCondEvent(Task,condition_in,phase_in)
%Must be modified if adding more tasks
%Can only use ONE phase as reference
switch upper(Task)%for now, keep it task specific
    case {'4POP'}
        numCond=2;
    case {'WM'}
        numCond=2;
end

%take the first numCond conditions only
condition_in=condition_in(1:numCond,:);
condition_in(:,1)=cellfun(@Number2Word,condition_in(:,1),...
    'UniformOutput',false);%make sure the names does not contain number
condition_struct=cell2struct(condition_in(:,2),condition_in(:,1),1);
for m = 1:size(condition_in,1)
    foobar=intersect(condition_struct.(...
        condition_in{m,1}),cell2mat(phase_in{2}));
    if ~isempty(foobar)
        cond{m}=foobar;
        cond_name{m}=condition_in{m,1};
    end
end
    

end

%%
function cond_ind=FindCurrentCondition(...
    epoch,urevent,list_to_look,look_dir)
%epoch is current epoch
%urevent is the entire list of EEG.urevent
%list_to_look is the list of potential things that current epoch may
%contain (e.g. conditions or accuracy information)
%look_direction:    backward (towards 1st event): 1 (default)
%                   forward (toward last event) : -1

%try use epoch.eventtype first
current_cond=intersect(cell2mat(epoch.eventtype),cell2mat(list_to_look));

if isempty(current_cond)%if event condition is not withint eventtype
    if nargin<4
        look_dir=1;%default lookup direction as backward
    end
    counter=1;%lookingup "counter" units back on the urevent array
    while 1
        %find index where it is going to look up
        IND=min(cell2mat(epoch.eventurevent))-(look_dir*counter);
        current_cond=intersect(urevent(IND).type,cell2mat(...
            list_to_look));
        if isempty(current_cond)
            counter=counter+1;
            continue;
        else
            break;
        end
    end
end

%now lookup which set of conditions that current epoch belongs to
foobar=cellfun(@(x) find(x==current_cond),list_to_look,...
    'UniformOutput',false);
cond_ind=find(~cellfun(@isempty,foobar));
end

%%
function ACC_code=FindACC(Accuracy,varargin)
%find target Accuracy code
if isempty(varargin)
    desired_acc='Correct';%default looking up Correct code
else
    desired_acc=varargin{1};
end
IND=strcmpi(Accuracy(:,1),desired_acc);
ACC_code=Accuracy{IND,2};
end


