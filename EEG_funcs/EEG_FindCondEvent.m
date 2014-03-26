function [cond,cond_name,numCond]=EEG_FindCondEvent(Task,condition_in,phase_in)
%Must be modified if adding more tasks
%Can only use ONE phase as reference
switch upper(Task)
    case {'4POP'}
        numCond=2;
    case {'WM'}
        numCond=3;
end

%take the first numCond conditions only
condition_in=condition_in(1:numCond,:);
condition_names=condition_in(:,1);%store condition names before conversion
condition_in(:,1)=cellfun(@Number2Word,condition_in(:,1),...
    'UniformOutput',false);
condition_struct=cell2struct(condition_in(:,2),condition_in(:,1),1);
for m = 1:size(condition_in,1)
    foobar=intersect(condition_struct.(...
        condition_in{m,1}),cell2mat(phase_in{2}));
    if ~isempty(foobar)
        cond{m}=foobar;
        cond_name{m}=condition_names{m,1};
    end
end
end