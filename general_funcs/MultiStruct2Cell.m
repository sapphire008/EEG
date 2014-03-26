function table_cell = MultiStruct2Cell(struct_array,table_form)
% table_cell = MultiStruct2Cell(struct_array,table_form)
% table_form: either long form (default) or short form
tmp = struct_array;%initiate first layer inspection
layer=1;%initiate first layer
array_leng = cell();
num_fields = cell();
field_list = cell();
while isstruct(tmp)
    [array_leng{layer},num_fields{layer},field_list{layer}] =...
        report_current_level(tmp);
    tmp_new = tmp(1).(field_list{layer}{1});
    clear tmp;
    tmp = tmp_new;
    clear tmp_new;
    layer = layer+1;
end

end

function [array_leng,num_fields,field_list] = ...
    report_current_layer(some_struct)
array_leng = length(some_struct);
field_list = fieldnames(some_struct);
num_fields = length(field_list);
end
 