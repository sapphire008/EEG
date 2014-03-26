function output_struct=CopyField(input_struct,field_names)
%pass input_struct.(field_names) to output_struct.(field_names)
%field_names are cell array of strings
%note that this will create a new structure
for n = 1:length(field_names)
    output_struct.(field_names{n})=input_struct.(field_names{n});
end
end