function [M1,M2,M3,M4]=absExtreme(V,varargin)
%[M1,M2,M3,M4] = absExtreme(V,...)
%M1: absolute max
%M2: index of absolute max
%M3: absolute min
%M4: index of absolute min
%
%[M1,M2] = absExtreme(V,...)
%if not specified to return max or min
%M1: [absolute_max,index_of_absolute_max]
%M2: [absolute_min,index_of_absolute_min]
%if specified to return max or min
%M1: absolute value
%M2: index of absolute value
%
%M = absExtreme(V,...)
%M: [M1,M2;M3,M4]


%find absolute extreme values of a vector
%Optional input:
%   'max': return absmax only
%   'min': return absmin only

if isempty(varargin)
    M=[absmax(V);absmin(V)];
else
    switch lower(varargin{1})
        case {'max'}
            M=absmax(V);
        case {'min'}
            M=absmin(V);
    end
end

%outputs
switch nargout
    case {0,1}
        M1 = M;
    case {2}
        if isempty(varargin)
            M1 = M(1,:);
            M2 = M(2,:);
        else
            M1 = M(1,1);
            M2 = M(1,2);
        end
case {4}
        if isempty(varargin)
            M1 = M(1,1);
            M2 = M(1,2);
            M3 = M(2,1);
            M4 = M(2,2);
        else
            error('Too many output requested');
        end
            
    otherwise
        error('Number of Output must be 1, 2, or 4');
end
end

function M_max=absmax(V)
%find absolute maximum
[~,I]=max(abs(V));
M_max=[V(I),I];
end

function M_min=absmin(V)
%find absolute minimum
[~,J]=min(abs(V));
M_min=[V(J),J];
end
