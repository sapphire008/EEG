function EEG_Out=EEG_rmbase(EEG,Task,Condition,phase_in,baseline_leng)
%Custom version of remove baseline. 
%Can remove baseline based on conditions of the epoch.

%find indices of epoch for each condition
[Vector_struct,~,~,cond_name]=EEG_makeVector(...
    EEG,Task,Condition,phase_in,'rmempty',1);

EEG_Out=EEG; %place holding

%remove baseline for each condition
for m = 1:length(cond_name)
    if isfield(Vector_struct.idx,cond_name{m})
        %select trials of current condition
        tmp_EEG=pop_select(EEG,'trial',Vector_struct.idx.(cond_name{m}));
        %remove baseline within current condition
        tmp_EEG=pop_rmbase(tmp_EEG,...
            [Condition{m,3}-baseline_leng Condition{m,3}]);
        %store the data after removing baseline back to EEG structure
        EEG_Out.data(:,:,Vector_struct.idx.(cond_name{m}))=...
            tmp_EEG.data;
        clear tmp_EEG;
    end
end
end






