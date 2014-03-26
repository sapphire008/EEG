function Comparison_Fields=EEG_makeComparisonFields(...
    Comparison_struct,level_list)
%clear;
%level_list = {{'SZ','C'},{'incong','cong','contrast'}};
% Comparison_table = {...
%     'comparison_type', 'Within',                      'Between';...
%     'Between_Groups',  {'incong','cong','contrast'},  {'SZ','C'};...
%     'Within_Groups',   {'SZ','C'},                    {'incong','cong'}};
% Comparison_struct = Table2Struct(Comparison_table,2);

comparison_types = fieldnames(Comparison_struct);
Comparison_Fields = struct();
for n =1:length(comparison_types)
    tmp = Comparison_struct.(comparison_types{n});
    Comparison_Fields.(comparison_types{n}) = cell(...
        length(tmp.Within),length(tmp.Between));
    switch sum(ismember(tmp.Within,level_list{1}))>0
        case {1}
            first = 1;
            second = 2;
            
        case {0}
            first = 2;
            second = 1;
    end
        
    for p = 1:length(tmp.Within)
        for q = 1:length(tmp.Between)
        Comparison_Fields.(comparison_types{n}){p,q}{first} = ...
            tmp.Within{p};
        Comparison_Fields.(comparison_types{n}){p,q}{second} = ...
            tmp.Between{q};
        end
    end
end
end