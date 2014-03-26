function str_out=ZeroPad(str_in,num,varargin)
% str_out=ZeroPad(str_in,num)
% str_in: the input string
% num: number of digits wanted to achieve
% Optional inputs
%   'PadPosition': 'start'(default)|'end'|[start index (double)]

pad_position=InspectVarargin(varargin,{'PadPosition'},{'start'});
zeros_str = repmat('0',1,num-length(str));

switch lower(pad_position)
    case {'start'}
        str_out = [zeros_str,str_in];
    case {'end'}
        str_out = [str,zeros_str];
    otherwise
        str_out = [str(1:(pad_position-1)),zero_str,str(pad_position:end)];    
end
end