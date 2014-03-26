function get_current=E4P_Conv2double(OBJ)
objClass=class(OBJ);%get the class of the object [char]
switch objClass
    case {'char'}%in case type 'str' or 'char'
        get_current=str2num(OBJ);
    case {'cell'}%in case of cell array
        if ischar(OBJ{1,1})
            get_current=str2num(OBJ{1,1});
        else
            get_current=OBJ{1,1};
        end
    otherwise %in case of double
        get_current=OBJ(1,1);
end
        