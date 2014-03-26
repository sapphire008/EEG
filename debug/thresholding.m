%thresholding
clear all;
subject = 'PFC102_032013';
task = 'WM';
load(['/nfs/jong_exp/EEG_PFC/subjects/',subject,'/',task,'/workspace.mat']);
eval(['!rm -r /nfs/jong_exp/EEG_PFC/subjects/',subject,'/',task,'/Data']);
eval(['!rm -r /nfs/jong_exp/EEG_PFC/subjects/',subject,'/',task,'/Spectrogram_gamma']);

for n = 1:size(Phase,1)
 EEG = pop_loadset('filepath',[params.Save_dir ,Phase{n,1} ...
        '/'],'filename',[params.subjectID,'-Analysis-Ready.set']);
 EEG=EEG_artifactRejection(EEG,2,params.artThresh);
    % Double check if the epoches are labeled with correct markers
    EEG=EEG_checkEpoch(EEG,'InvalidMarker',255,'DisplayEpoch',0);
    %save after final artifact reject and marker check
    EEG = pop_saveset(EEG,'filepath',[params.Save_dir ,Phase{n,1} ...
        '/'],'filename',[params.subjectID,'-Analysis-Ready.set']);
    numEpoch.(Phase{n,1})=length(EEG.epoch);
    
    save([params.Save_dir,'workspace.mat'],'elecs','bad_comps',...
    'numEpoch','-append');
end