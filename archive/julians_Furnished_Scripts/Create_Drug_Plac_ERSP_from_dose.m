%% change dose1/dose2
% make sure there are 2 ERSPs, named 'ERSP_group_doseX'
%manually save after. Don't forget freqs, timesout

%set group (mc or ms)
group = 'ms'
dose = 'dose3' %this only matters for the dose3 ms

if ~isequal(dose,'dose3') %for dose1 and dose2

if isequal(group,'mc')
    dose1_is_Drug = [2 4 5 6 7 14 16 17 18 20 23 103 110 112] %mc
end

if isequal(group,'ms')    
    dose1_is_Drug = [4 6 10 12 17 18 20 23 24 29 30 32 33 36] %ms
end

%get ERSP size
eval(['max_subj = length(ERSP_',group,'_dose1.dose1)'])

for s = [1:max_subj]
   
       if sum(dose1_is_Drug==s)==1 %if dose1 was Drug day
        
       eval(['ERSP.Drug(s) = ERSP_',group,'_dose1.dose1(s)'])
       eval(['ERSP.Plac(s) = ERSP_',group,'_dose2.dose2(s)'])     
           
       else %if dose1 was Plac day
        eval(['ERSP.Drug(s) = ERSP_',group,'_dose2.dose2(s)'])
        eval(['ERSP.Plac(s) = ERSP_',group,'_dose1.dose1(s)'])     
                
       
    end
end

else
    
    dose3_is_Drug = [4 7 10 14 16 17 20 23 24 27 28 29 33 35 36 37 30]
    dose3_is_Plac = [2 8 12 13 15 18 22 25 32 38]
    
    ERSP.Drug(dose3_is_Drug) = ERSP_ms_dose3.dose3(dose3_is_Drug)
    ERSP.Plac(dose3_is_Plac) = ERSP_ms_dose3.dose3(dose3_is_Plac)
end
    
    
    




