% set the paths and variables so they match the group

if s >= 100
    zeroes = ''
else

    if s >= 10
        zeroes = '0'
      

    else

        zeroes = '00'
       
    end

end


cd(preprocdir)
cd([num2str(group),num2str(zeroes),num2str(s)])
EEG = pop_loadset( 'filename', [num2str(group),num2str(zeroes),num2str(s),'-',dose,'-PostR.set']);

EEG = pop_runica( EEG, 'icatype', 'runica','extended',1,'pca',75);

EEG = pop_saveset(EEG, 'filename', [num2str(group),num2str(zeroes),num2str(s),'-',dose,'-PostR-ICA.set']);
clear EEG

eeglab

