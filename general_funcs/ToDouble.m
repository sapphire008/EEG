function get_current=ToDouble(OBJ)
%converting an EEG.epoch.eventtype to a double
objClass=class(OBJ);%get the class of the object [char]
switch objClass
    case {'char'}%in case type 'str' or 'char'
        get_current=str2num(OBJ);
    case {'float','uint16','uint8','double','logical'}%numeric
        get_current=double(OBJ);
    case {'cell'}%in case of cell array
        if ischar(OBJ{1,1})%char inside a cell
            get_current=cellfun(@str2num,OBJ);
        else
            get_current=cell2mat(OBJ);
        end
    
end
        