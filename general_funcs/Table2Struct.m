function table_struct=Table2Struct(table_cell,dim)
%table_struct = TABLE2STRUCT(table_cell,dim)
%converting a NxM table of cell array in to a structure with N or M fields,
%and within each field has M or N subfields, respectively.
%dim specifies which dimension will be the primary field
%       F1    F2    F3    F4
% f1    V1    V2 ...
% f2    .     .
% f3    .     .
% f4    .     .
%
% <-=-> equivalently,
%
% table_cell.f2.F1=V5; %when dim == 1
% OR
% table_cell.F1.f1=V1; %when dim == 2

[numRow,numCol]=size(table_cell);

for i=2:numRow
    for j = 2:numCol
        switch dim
            
            case {1}
                table_struct.(table_cell{1,j}).(table_cell{i,1})=...
                    table_cell{i,j};
                
            case {2}
                table_struct.(table_cell{i,1}).(table_cell{1,j})=...
                    table_cell{i,j};
                
        end
    end
end

end