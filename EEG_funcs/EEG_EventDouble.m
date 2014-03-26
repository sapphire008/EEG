function get_current=EEG_EventDouble(OBJ)
%converting an EEG.epoch.eventtype to a double
objClass=class(OBJ);%get the class of the object [char]
switch objClass
    case {'char'}%in case type 'str' or 'char'
        get_current=str2num(OBJ);
    case {'float','uint16','uint8','double','logical'}%numeric
        get_current=double(OBJ);
    case {'cell'}%in case of cell array
        if ischar(OBJ{1,1})
            get_current=str2num(OBJ{1,1});
        else
            get_current=OBJ{1,1};
        end
    
end
        