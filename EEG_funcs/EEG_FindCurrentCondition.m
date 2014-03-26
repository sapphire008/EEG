function cond_ind=EEG_FindCurrentCondition(...
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