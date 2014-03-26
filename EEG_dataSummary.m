%get data quanlity
base_dir = '/nfs/jong_exp/EEG_PFC/subjects/';
save_dir = '/nfs/jong_exp/EEG_PFC/GroupAnalysis/';
subjects = {'PFC101_031813','PFC102_032013','PFC103_040313',...
    'PFC104_041213','PFC200_030813','PFC201_032813','PFC202_032213',...
    'PFC203_040813','PFC204_041713'};
tasks = {'4POP','WM'};
tmp_tasks = {'FourPOP','WM'};
Phase = {'Cue','Probe','Whole'};
Conditions.FourPOP = {'Green','Red'};
Conditions.WM = {'TwoFace','OneFace'};
Data_Summary=cell2struct(cell(1,length(tasks)),tmp_tasks,2);

for tk = 1:length(tasks)
    for s = 1:length(subjects)
        if ~exist([base_dir,subjects{s},'/',tasks{tk}],'dir')
            continue;
        else
            %store subject name of current taks
            Data_Summary.(tmp_tasks{tk})(s).subject = subjects{s};
        end
        for n = 1:length(Phase)
            for c = 1:length(Conditions.(tmp_tasks{tk}))
                clear EEG;
                EEG = pop_loadset([base_dir,subjects{s},'/',tasks{tk},...
                    '/',Phase{n},'/Conditions/',subjects{s},'-',...
                    Phase{n},'_',Conditions.(tmp_tasks{tk}){c},'_all.set']);
                
                %summarize data
                Data_Summary.(tmp_tasks{tk})(s).(Phase{n}).(...
                    Conditions.(tmp_tasks{tk}){c}).num_elecs = EEG.nbchan;
                Data_Summary.(tmp_tasks{tk})(s).(Phase{n}).(...
                    Conditions.(tmp_tasks{tk}){c}).num_comps = size(EEG.icaweights,1);
                Data_Summary.(tmp_tasks{tk})(s).(Phase{n}).(...
                    Conditions.(tmp_tasks{tk}){c}).num_trials = length(EEG.epoch);
            end
        end
    end
end
save([save_dir,'DataSummary.mat'],'Data_Summary');