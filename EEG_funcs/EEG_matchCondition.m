function RawDataSet=EEG_matchCondition(EEG,Phase,Condition,RawDataSet)
%under the condition of 4POP, i.e. red corersponds to incongreuent, and
%green corresponds to congruent
%find which condition that each trial/epoch/block belongs to

for ep = 1:1:length(EEG.epoch)
     %current_condition is the event code formatted as double
     current_condition=get_current_condition(EEG,ep,Phase{n,1});
    %storing what epoch #s belong to this condition
    switch current_condition%PROBLEM!!
        case {5}%green
        RawDataSet.(Phase{n,1}).idx.(Condition{1,1})(end+1) = ep; 
        case {6}%red
        RawDataSet.(Phase{n,1}).idx.(Condition{2,1})(end+1) = ep;
    end
end


end

%get the current condition given current phase
function current_condition=get_current_condition(EEG,ep,Condition)
cond_list=cell2mat(Condition(2,:)');
current_event=ToDouble(EEG.epoch(ep.eventtype));
%to see if current event codes has condition information
event_in_cond=ismember(cond_list,current_event);
switch sum(event_in_cond)>1E-6
    case {0}%if condition does not occur in the event type
        urevent_IND=ToDouble(EEG.epoch(ep).eventurevent{1,1})-1;
             current_condition=EEG_EventDouble(...
                EEG.urevent(urevent_IND).type);     
    case {1}%if condition occurs in event type
        current_condition=cond_list(event_in_cond>1E-6);
        
end
end