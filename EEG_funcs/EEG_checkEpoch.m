function EEG=EEG_checkEpoch(EEG,varargin)
% This function keeps all the events markers as numeric
%default invalid marker: 255
flag=InspectVarargin(varargin,...
    {'InvalidMarker','DisplayEpoch'},...
    {[255,NaN],0});
disp(['Performing additional epoch check to make sure',...
    ' no invalid or duplicated markers inserted']);
% correct epoch field
for m = 1:length(EEG.epoch)
    %check if the eventtype are numeric, if not need to convert everything
    %to numeric and find out what is causing the eventtype to be
    %non-numeric
    EEG.epoch(m)=removeNonNumeric(EEG.epoch(m));
    %find index of invalid markers
    InvalidNumMarkerInd=...
        checkInvalidNumMarker(EEG.epoch(m),flag.InvalidMarker);
    %find index of duplicate markers
    DuplicateMarkerInd=checkDuplicateMarker(EEG.epoch(m));
    %find index of marker as NaN
    NaNMarkerInd = chekcNaN(EEG.epoch(m));
    marker_keeping_vect=setdiff(1:length(EEG.epoch(m).event),...
        union(union(InvalidNumMarkerInd,DuplicateMarkerInd),NaNMarkerInd));
    EEG.epoch(m)=structfun(@(x) x(marker_keeping_vect),EEG.epoch(m),...
        'UniformOutput',false);
    if flag.DisplayEpoch
        disp([m,':',EEG.epoch(m).eventtype]);
    end
end

%correct urevent and event fields
EEG.urevent = correctNonDoubleEvent(EEG.urevent,'type');
EEG.event = correctNonDoubleEvent(EEG.event,'type');

end


function current_epoch=removeNonNumeric(current_epoch)
%check if all is numeric
tmp=~cellfun(@isnumeric,current_epoch.eventtype);
if sum(tmp)<1%this means all the eventtype are numeric
    return;
else
    %check if they are numeric strings
    tmp2=cellfun(@(x) isempty(regexp(x,'\D')),...
        current_epoch.eventtype,'UniformOutput',false); %#ok<RGXP1>
    %remove markers that are non-numeric char
    current_epoch=structfun(@(x) x(cell2mat(tmp2)>0),current_epoch,...
        'UniformOutput',false);
    %convert numeric char markers into numerics
    current_epoch.eventtype=cellfun(@(x) str2double(x),...
        current_epoch.eventtype,'UniformOutput',false);
    %also need to correct urevent and event field of the dataset
end
end

function Event=correctNonDoubleEvent(Event,field_str)
%convert a field of numeric char to a double
for k = 1:length(Event)
    if ischar(Event(k).(field_str))
        %will convert numeric char to the corersponding number and
        %non-numeric chars into NaN, which is okay
        Event(k).(field_str)=str2double(Event(k).(field_str));
    end
end
end

function  InvalidNumMarkerInd=...
    checkInvalidNumMarker(current_epoch,InvalidMarker)
%find index of invalid numerical markers
[~,InvalidNumMarkerInd,~]=...
    intersect(cell2mat(current_epoch.eventtype),InvalidMarker);
end

function DuplicateMarkerInd=checkDuplicateMarker(current_epoch)
%each phase should have only one type of marker, unless this rule will be
%modified later, then specify the number of markers should be present
%at each phase
et_vect=cell2mat(current_epoch.eventtype);
tmp_unique=unique(et_vect);
    if length(tmp_unique)~=length(et_vect)
        N=hist(et_vect,tmp_unique);
        tmp_dup_ind=find(et_vect==tmp_unique(N>1));
        DuplicateMarkerInd=tmp_dup_ind(2:end);%keep ind>1 as duplicate ind
    else
        DuplicateMarkerInd=[];
    end
end

function NaNInd = chekcNaN(current_epoch)
NaNInd = find(isnan(cell2mat(current_epoch.eventtype)));
end
