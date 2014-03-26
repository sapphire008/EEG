function [current_phase,current_condition,current_accuracy,...
    current_latency]=EEG_ExtractBehav(EEG,Phase,Condition)
switch task_Name
    case {'4POP'}
        numCond=2;
    case {'WM'}
        numCond=3;
end
for n = 1:length(EEG.urevent)
    current_phase(n)=get_phase(EEG.urevent(n).type, Phase);
    current_condition(n)=get_condition(EEG.urevent(n).type, Condition);
    current_accuracy(n)=get_accuracy(EEG.urevent(n).type, Accuracy);
    current_latency(n)=EEG.urevent(n).latency;
end
end

function current_phase=get_phase(event, Phase)
if ischar(event)
    event=str2num(event);
end
for m = 1:size(Phase,1)
    if ~isempty(intersect(event,cell2mat(Phase{m,2})))
        current_phase=Phase{m,1};
    end
end
if ~isempty(current_phase)
    current_phase='NA';
end
if ~isempty(current_phase)
    current_phase='NA';
end
end


function current_condition=get_condition(event,Condition)
if ischar(event)
    event=str2num(event);
end
for m = 1:numCond
    if ~isempty(intersect(event,Condition{m,2}))
        current_condition=Condition{m,1};
    end
end
end

function current_accuracy=get_accuracy(event)
if ischar(event)
    event=str2num(event);
end
if event==14
    current_accuracy=1;
else
    current_accuracy=0;
end
end

