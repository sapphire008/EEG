%% reject ICs

if s >= 100
        zeroes = ''
    else

        if s >= 10
            zeroes = '0'

        else

            zeroes = '00'
        end

    end

    cd(preprocdir);   
    
EEG = pop_loadset( 'filename', [num2str(group),num2str(zeroes),num2str(s),'-',dose,'-PostR-ICA.set'], 'filepath', [preprocdir,num2str(group),num2str(zeroes),num2str(s)]);

load(IC_file)
eval(['comps = ',group,'_removed_components_dose']);

load(YG_file) 
 eval(['comps_YG = ',group,'_YG_microsaccades_dose']);

if isequal(dose,'dose1') || isequal(dose,'dose-1') || isequal(dose,'dose-1')
    reject = unique([comps{s,2} comps_YG{s,2}])       
end
if isequal(dose,'dose2')
    reject = unique([comps{s,3} comps_YG{s,3}])
elseif isequal(dose,'dose3') 
    reject = unique([comps{s,4} comps_YG{s,4}])
end


EEG = pop_subcomp(EEG,reject,0);

EEG = pop_saveset(EEG, 'filename', [num2str(group),num2str(zeroes),num2str(s),'-',dose,'-Analysis2.set'],'filepath', [preprocdir,num2str(group),num2str(zeroes),num2str(s)]);

clear EEG

EEG=pop_loadset('filename',[num2str(group),num2str(zeroes),num2str(s),'-',dose,'-Analysis2.set'], 'filepath', [preprocdir,num2str(group),num2str(zeroes),num2str(s)]);

%EEG = pop_select( EEG, 'nochannel', {'62' '82'});

EEG = pop_eegthresh(EEG,1,[1:EEG.nbchan],-50,50,EEG.xmin,EEG.xmax,1,0);


    %Reject marked epochs and save the dataset
EEG = eeg_rejsuperpose( EEG, 1, 1, 1, 0, 0, 0, 0, 0);
EEG = pop_rejepoch( EEG, find(EEG.reject.rejglobal), 0);
EEG = pop_saveset(EEG, 'filename', [num2str(group),num2str(zeroes),num2str(s),'-',dose,'-Analysis2-Ready.set'],'filepath', [preprocdir,num2str(group),num2str(zeroes),num2str(s)]);

clear EEG
eeglab redraw


