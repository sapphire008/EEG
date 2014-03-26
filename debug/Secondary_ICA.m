%secondary runica
base_dir='/nfs/jong_exp/EEG_PFC/subjects/PFC200_030813/4POP/';
for n = 1:size(Phase,1)
    clear EEG;
    %EEG_ICA = pop_loadset([base_dir,Phase{n,1},'/PFC200_030813-Post-ICA.set']);
    EEG = pop_loadset([base_dir,Phase{n,1},'/PFC200_030813-epoched-PostR-visual-inspection.set']);
    %EEG.data = [];
    %EEG.data = EEG_ICA.data;
    %clear EEG_ICA;
    EEG = pop_runica(EEG, 'icatype', 'runica','extended',1);
    EEG = pop_saveset(EEG,...
        'filepath',[params.Save_dir Phase{n,1} '/'],...
        'filename',[params.subjectID,'-ICA-ready.set']);
end
%%
clear n;
for n = 1:length(Phase)
    EEG = pop_loadset([base_dir,Phase{n,1},'/PFC201_032813-secondary-ICA-Ready.set']);
    [EEG,bad_comps.(Phase{n,1})]=EEG_Rejector(EEG,2,'numComp',20,...
        'PlotOff',0);
    %reject the ICA components and save
    EEG = pop_saveset(EEG,...
        'filepath',[params.Save_dir Phase{n,1} '/'],...
        'filename',[params.subjectID,'-secondary-Post-ICA.set']);
    
    
end