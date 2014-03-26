function cell_array=FillEmptyCell(cell_array,varargin)
% Fill empty cell of a cell array with a matrix of numbers
% default fill with nan matrix of the same size as its previous
% neighbor.
% Optional Inputs:
%       'MatrixDimension': vector of matrix dimension, e.g. [5 3] = 5 x 3
%       'FillValue': value to fill, input as string, default 'NaN'
flag=InspectVarargin(varargin,{'MatrixDimension','FillValue',},...
    {'neighbor','NaN'});

%find empty cell index

[emptyIND,refIND]=findEmptyIND(cell_array);

for m = 1:length(emptyIND);
    %find the dimension of fill matrix
    switch flag.MatrixDimension
        case {'neighbor'}%get neighbor matrix size
            mat_dim=size(cell_array{refIND(m)});
        otherwise
            mat_dim=flag.MatrixDimension;
    end
    
    %fill each empty cell with the matrix specified
    switch flag.FillValue
        case {'0'}
            cell_array{emptyIND(m)}=zeros(mat_dim);
        case {'NaN'}
            cell_array{emptyIND(m)}=NaN(mat_dim);
        otherwise
            cell_array{emptyIND(m)}=flag.FillValue*ones(mat_dim);
    end
    
    %check data type consistency only in neighbor mode
    if strcmpi(flag.MatrixDimension,'neighbor')
            datatype=class(cell_array{refIND(m)});
            cell_array{emptyIND(m)}=reshape(cast(cell_array{...
                emptyIND(m)}(:),datatype),mat_dim);
    end       
end
end


function [emptyIND,refIND]=findEmptyIND(cell_array)
%find indices of nearest non-empty cells of empty cells
%in case there are two nearest neighbors (forward and backward), use
%backward nearest neighbor
emptyIND=find(cellfun(@isempty,cell_array));
nonempty_IND=find(~cellfun(@isempty,cell_array));
refIND = zeros(1,length(emptyIND));
for m = 1:length(emptyIND)
    [~,tmp_IND] = absExtreme(nonempty_IND-emptyIND(m),'min');
    refIND(m) = nonempty_IND(tmp_IND);
end
end





















