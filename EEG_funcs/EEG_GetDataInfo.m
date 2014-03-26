function DataInfo=EEG_GetDataInfo(...
    base_dir,Task_list,tmpl_subj)
%this can only be done after the processing paradigms are determined and at
%least one subject is done

% base_dir='/nfs/jong_exp/EEG_PFC/';
% Task_list={'4POP','WM'};
% tmpl_subj='PFC103_040313';
%place holding task information
DataInfo=struct();
for tk = 1:length(Task_list)
    %use one subject as a template to get params
    tmpl_name=getfield(dir([base_dir,'subjects/',tmpl_subj,'/',...
        Task_list{tk},'/Data/ERSP_PFC*.mat']),'name');
    %load ERSP <Phase x Condition x Elecs x Frequency x Time> structure
    load([base_dir,'subjects/',tmpl_subj,'/',Task_list{tk},...
        '/Data/',tmpl_name]);
    load([base_dir,'subjects/',tmpl_subj,'/',Task_list{tk},...
        '/workspace.mat'],'TFparams');
    
    %get structure info
    Phases=fieldnames(ERSP);
    Conditions=fieldnames(ERSP.(Phases{1}));
    Elecs=length(ERSP.(Phases{1}).(Conditions{1}).lead);
    
    %passing data to the DataInfo structure and return the structure
    Task_list{tk}=Number2Word(Task_list{tk});
    DataInfo.(Task_list{tk}).Phases=Phases;
    DataInfo.(Task_list{tk}).Conditions=Conditions;
    DataInfo.(Task_list{tk}).Elecs=Elecs;
    
    %get freqs x time info
    for ph = 1:length(Phases)
        [DataInfo.(Task_list{tk}).(Phases{ph}).Freqs,...
            DataInfo.(Task_list{tk}).(Phases{ph}).Times]=...
            size(ERSP.(Phases{ph}).(Conditions{1}).lead{1});
        DataInfo.(Task_list{tk}).(Phases{ph}).times_out = TFparams.(Phases{ph}).times_out;
        DataInfo.(Task_list{tk}).(Phases{ph}).freqs_out = TFparams.(Phases{ph}).freqs_out;
    end
end
end