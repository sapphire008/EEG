function [S,empty_ind]=RemoveEmptyField(S)
field_name_list=fieldnames(S);%get a list of field names
empty_ind=structfun(@(x) isempty(x),S);%get the index of empty field
S=rmfield(S,field_name_list(empty_ind));%remove empty field
end