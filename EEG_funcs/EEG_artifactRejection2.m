function E4P_artifactRejection2(preprocdir,subjectID, phase,rejType)
%rejType: 1.Primary artifact reject: based on probability
%         2.Primary artifact reject: based on voltage deflection
%         3.Final artifact reject

loadName='';%place holding name of the file to be loaded
saveName='';%place holding name of the file to be saved
numLoop=1;%default no loop in primary reject

% Set parameters for each case of rejType
switch rejType
    case {1,2}
        loadName='-epoched.set';
        saveName='-epoched-PostR.set';
        thresh=150;%used in reject by voltage deflection
        disp('case 1 or 2 is used');%for debug only
    case {3}
        loadName='-Post-ICA.set';
        saveName='-Analysis-Ready.set';
        thresh=50;%used in reject by voltage deflection
        numLoop=length(phase(:,1));%is number of phases in final reject
        disp('case 3 is used');%for debug only
end%end of switch

% Executing artifact reject
for i=1:numLoop
    EEG = pop_loadset('filename',[subjectID,loadName], ...
        'filepath',[preprocdir, subjectID, '/', phase{i,1}, '/']);
    switch rejType
        case {1}
            % Reject epochs based on probability
            [EEG, locthresh, globthresh, nrej] = pop_jointprob(EEG,...
                1,[1:EEG.nbchan] ,3,5,1,0,0);
            %Reject marked epochs and save the dataset
            EEG = eeg_rejsuperpose( EEG, 1, 1, 1, 1, 1, 1, 1, 1);
            EEG = pop_rejepoch( EEG, find(EEG.reject.rejglobal), 0);
            disp('case 1 is used');%for debug only
        case {2,3}
            % Reject based on voltage deflection
            EEG = pop_eegthresh(EEG,1,[1:EEG.nbchan],-thresh,thresh,...
                EEG.xmin,EEG.xmax, 0);
            disp('case 2 or 3 is used');%for debug only
    end%end of switch
    % Save dataset
    EEG = pop_saveset( EEG,'filename',[subjectID,saveName],...
        'filepath',[preprocdir, subjectID, '/', phase{i,1}, '/']);
    
    clear EEG locthresh globthresh nrej
end


% Oiriginal code separated by case
% switch rejType
%     case {1}
%         % Primary artifact reject 
%         EEG = pop_loadset([subjectID,'-epoched.set'], [preprocdir, subjectID, '/', phase{i,1}, '/'])
%         
%         % Reject epochs based on probability (default)
%         [EEG, locthresh, globthresh, nrej] = pop_jointprob(EEG,1,[1:EEG.nbchan] ,3,5,1,0,0);
%         %Reject marked epochs and save the dataset
%         EEG = eeg_rejsuperpose( EEG, 1, 1, 1, 1, 1, 1, 1, 1);
%         EEG = pop_rejepoch( EEG, find(EEG.reject.rejglobal), 0);
%         EEG = pop_saveset( EEG,[subjectID,'-epoched-PostR.set'], [preprocdir, subjectID, '/', phase{i,1}, '/'])
%     case {2}
%         % Reject based on voltage deflection
%         EEG = pop_loadset([subjectID,'-epoched.set'], [preprocdir, subjectID, '/', phase{i,1}, '/'])
%         EEG = pop_eegthresh(EEG,1,[1:EEG.nbchan],-150,150,EEG.xmin,EEG.xmax,1,0);
%         % save dataset
%         cd([preprocdir,subjectID])
%         EEG = pop_saveset( EEG,[subjectID,'-epoched-PostR.set'], [preprocdir, subjectID, '/', phase{i,1}, '/'])
%     case {3}
%         % Final Artifact Reject
%         for i = 1:numPhase
%             EEG = pop_loadset( [subjectID,'-Post-ICA.set'], [preprocdir, subjectID, '/', phase{i,1}, '/'])
% 
%             EEG = pop_eegthresh(EEG,1,[1:EEG.nbchan],-50,50,EEG.xmin,EEG.xmax,1,0);
% 
%             EEG = pop_saveset( EEG,'filename',[subjectID,'-Analysis-Ready.set'], 'filepath',[preprocdir, subjectID, '/', phase{i,1}, '/'])
%         end %end of for-loop
% end %end of switch case

end%end of function
