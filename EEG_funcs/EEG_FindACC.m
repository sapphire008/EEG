function ACC_code=EEG_FindACC(Accuracy,varargin)
%find target Accuracy code
if isempty(varargin)
    desired_acc='Correct';%default looking up Correct code
else
    desired_acc=varargin{1};
end
IND=strcmpi(Accuracy(:,1),desired_acc);
ACC_code=Accuracy{IND,2};
end