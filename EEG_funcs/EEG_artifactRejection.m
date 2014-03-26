function [EEG,IND] = EEG_artifactRejection(EEG, rejType,varargin)
%EEG = E4P_artifactRejection2(EEG,rejType,thresh)
% rejType:
    %1: reject by probability
    %2: reject by voltage deflection
% thresh: necessary only when using type 2.
flag = InspectVarargin(varargin,{'thresh','rejEpoch'},{50,1});

switch rejType%type double
    case {1}%reject by probability
        [EEG, locthresh, globthresh, nrej] = pop_jointprob(EEG,...
                1,[1:EEG.nbchan] ,3,5,1,0,0);
        IND = 0;
        %Reject marked epochs and save the dataset
        EEG = eeg_rejsuperpose( EEG, 1, 1, 1, 1, 1, 1, 1, 1);
        EEG = pop_rejepoch( EEG, find(EEG.reject.rejglobal), 0);
        disp('case 1:reject by probability is used');%for debug
    case {2}%reject by voltage deflection
            [EEG,IND] = pop_eegthresh(EEG,1,[1:EEG.nbchan],...
                -flag.thresh,flag.thresh,EEG.xmin,EEG.xmax,1,...
                flag.rejEpoch);
            disp('case 2:reject by voltage deflection is used');%for debug

end      
end