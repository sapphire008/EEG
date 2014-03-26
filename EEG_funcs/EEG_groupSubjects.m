function subj_vect = EEG_groupSubjects(subj_list,subj_group,study_name)
%get the numbering of the subject and separate them into group, that is
%within the number_range
subj_vect = struct();
max_num_len=max(cellfun(@(x) max(floor(log10(x))+1), subj_group(:,2)));
tmp_str=cellfun(@(x) x(length(study_name)+(1:max_num_len)),...
    subj_list,'UniformOutput',false);
subj_num=cell2mat(cellfun(@(x) str2double(x),tmp_str,'UniformOutput',false));

for n = 1:size(subj_group,1)
    subj_vect.(subj_group{n,1})=find(subj_num<=max(subj_group{n,2}) & ...
        subj_num>=min(subj_group{n,2}));
end
end