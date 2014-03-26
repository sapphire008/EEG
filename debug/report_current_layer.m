function [array_leng,num_fields,field_list] = ...
    report_current_layer(some_struct)
array_leng = length(some_struct);
field_list = fieldnames(some_struct);
num_fields = length(field_list);
end