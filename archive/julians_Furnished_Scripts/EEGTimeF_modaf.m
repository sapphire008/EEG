%%%% wavelet the data using newtimef
%Author: Glenn Gomes, Julian Cheng
%
%Changelog:
%   1/22/2013:  Found a bug that was causing excessive NaNs in the result,
%               found to be caused by the pre-wavelet algorithm not
%               recognizing the data type stored in the EEG structure. This
%               issue is fixed and code added to warn if either cong_idx or
%               incong_idx is empty.
%   1/25/2013:  Added code to allow disabling of baseline subtraction, used
%               to diagnose if baseline subtraction is working 
%   1/28/2013:  Added code to allow switching between different baseline
%               subtraction methods
%   2/4/2013:   Added code to allow enabling logrithmic scaling of
%               frequency
%   2/5/2013:   Added code to allow enabling of permutation
%   3/19/2013:  Added code to allow for block-wise trial sorting, instead
%               of just all green vs all red. Note: this algorithm assumes
%               that each block has 82 trials and there are 8 blocks.
%   3/28/2013:  Added recording of History and source-id in output file,
%               which is adapted back from the SC scripts. Note that all
%               files processed after this date will have this included by
%               default. Also added file search for "dose1" files being
%               named to "dose-1". Also added a check to ensure output
%               directory exists, as well as overwrite confirmation

%default debug settings
if ~exist('DEBUG_UseDS1ForMCPC','var')
    DEBUG_UseDS1ForMCPC = false;
end
if ~exist('DEBUG_DisableBaseline','var')
    DEBUG_DisableBaseline = false;
end
if ~exist('DEBUG_FLAG_BaselineMethod','var')
    DEBUG_FLAG_BaselineMethod = 0;
end
if ~exist('DEBUG_UseLogrithmicFrequencyScaling','var')
    DEBUG_UseLogrithmicFrequencyScaling = false;
end
if ~exist('DEBUG_UsePermutation','var')
    DEBUG_UsePermutation = false;
end
if ~exist('DEBUG_intAlpha','var')
    DEBUG_intAlpha = 0.05;
end

%flag definitions
if ~exist('boolUseBlockSeperation','var')
    boolUseBlockSeperation = false;
    lstBlocks2Process = [];
end

%construct history
History = fnGetHistory('id',fnGetUniqueString());

%check if output directory exists
if ~exist(dirOutput,'dir')
    fprintf('Warning: output folder does not exist and will be automatically created\n\n')
    [boolSuccess,strMessage,strMessageID] = mkdir(dirOutput);
    if ~boolSuccess
        error('%s: %s',strMessageID,strMessage)
    end
end
clear boolSuccess strMessage strMessageID

%check if destination file already exists
if exist([dirOutput,'/ERSP_',num2str(ends),'.mat'],'file') || exist([dirOutput,'/ITC_',num2str(ends),'.mat'],'file')
    fprintf('Warning: target file already exists\n\n');
    beep;
    strUserInput = input('Overwrite existing file? (Y/N): ','s');
    if isempty(regexpi(strUserInput,'y'))
        error('Operation aborted')
    end
end
clear strUserInput

lstSubjectsNotFound = [];
for s = [sbjs]
    if s >= 100
        zeroes = '';
    else

        if s >= 10
            zeroes = '0';

        else

            zeroes = '00';
        end
    end
    
    cd(preprocdir);
    
    %update 3/28/2013
    %checks if subject folder exists
    if exist([num2str(group),num2str(zeroes),num2str(s)],'dir')
        cd([num2str(group),num2str(zeroes),num2str(s)])
    else
        lstSubjectsNotFound(end+1) = s;
        continue
    end
        
    %update 3/28/2013
    %checks if the target file exists, otherwise inserts a dash between the
    %alpha and numeric portions of "dose"
    strAlpha = dose(isstrprop(dose,'alpha'));
    strNumeric = dose(isstrprop(dose,'digit'));
    if exist([num2str(group),num2str(zeroes),num2str(s),'-',dose,'-Analysis2-Ready.set'],'file')
        EEG=pop_loadset('filename',[num2str(group),num2str(zeroes),num2str(s),'-',dose,'-Analysis2-Ready.set'],'filepath', [preprocdir,num2str(group),num2str(zeroes),num2str(s)]);
    elseif exist([num2str(group),num2str(zeroes),num2str(s),'-',strAlpha,'-',strNumeric,'-Analysis2-Ready.set'],'file')
        EEG=pop_loadset('filename',[num2str(group),num2str(zeroes),num2str(s),'-',strAlpha,'-',strNumeric,'-Analysis2-Ready.set'],'filepath', [preprocdir,num2str(group),num2str(zeroes),num2str(s)]);
    else
        lstSubjectsNotFound(end+1) = s;
        continue
    end
    clear strAlpha strNumeric
        
    elecs = {EEG.chanlocs.labels};

%% Cue

if isequal(phase, 'Cue')

    if isequal(group,'mc')
    
        if (s <= 7) && DEBUG_UseDS1ForMCPC %for DS1
            event_num_incong = 10;
            event_str_incong = '10';
            event_num_cong = 20;
            event_str_cong = '20';
        else
            event_num_incong = 11;
            event_str_incong = '11';
            event_num_cong = 21;
            event_str_cong = '21';
        end

    else 
        event_num_incong = 11;
        event_str_incong = '11';
        event_num_cong = 21;
        event_str_cong = '21';
    end

    cong_idx = [];
    incong_idx = [];

    for epoch_num = 1:length(EEG.epoch)

        if ischar(EEG.epoch(1,epoch_num).eventtype)
            if EEG.epoch(1,epoch_num).eventtype == event_str_incong
                incong_idx = [incong_idx; epoch_num];
            end
            if EEG.epoch(1,epoch_num).eventtype == event_str_cong
                cong_idx = [cong_idx; epoch_num];
            end
        end
        
        if iscell(EEG.epoch(1,epoch_num).eventtype)
            if ischar(EEG.epoch(epoch_num).eventtype{1})
                if EEG.epoch(epoch_num).eventtype{1} == event_str_incong
                    incong_idx = [incong_idx; epoch_num];
                end
                if EEG.epoch(epoch_num).eventtype{1} == event_str_cong 
                    cong_idx = [cong_idx; epoch_num];
                end
            else
                if EEG.epoch(epoch_num).eventtype{1} == event_num_incong  % normally, you'd use this line.
                    incong_idx = [incong_idx; epoch_num];
                end
                if EEG.epoch(epoch_num).eventtype{1} == event_num_cong 
                    cong_idx = [cong_idx; epoch_num];
                end
            end
        end
                
        if isnumeric(EEG.epoch(1,epoch_num).eventtype)   
            if EEG.epoch(epoch_num).eventtype == event_num_incong  % normally, you'd use this line.
                incong_idx = [incong_idx; epoch_num];
            end
            if EEG.epoch(epoch_num).eventtype == event_num_cong 
                cong_idx = [cong_idx; epoch_num];
            end
        end
    end
    
    %Update: 1/22/2013
    %data integrity check
    if isempty(incong_idx)
        error('Failed to find first condition (incongruent) for subject %i cue',s)
    elseif isempty(cong_idx)
        error('Failed to find second condition (congruent) for subject %i cue',s)
    end
end

%% Probe

if isequal(phase, 'Probe')

    cong_idx = [];
    incong_idx = [];
    none_idx = [];

    for i = 1:length(EEG.epoch)

        if isequal(group,'mc')

            if s<8 %%for DS1
                trial_number = EEG.epoch(1,i).eventurevent(1);
                event_code = EEG.epoch(1,i).eventtype(1);  %not sure why this is used (JC)
            end

            if s>7
                %Update: 1/22/2013
                %added check to avoid problems with numeric type event codes
                %Note: this logic is undocumented, not sure why trial
                %number adjustment is necessary
                %Update: 3/19/2013
                %Note: found what this logic is for. It seems that for
                %later subjects the urevent doesn't match up with the
                %behave count anymore; somehow Glenn figured out a formula
                %to fix this issue and get the real urevent number back
                if isnumeric(EEG.epoch(1,i).eventurevent)
                    trial_number = (EEG.epoch(1,i).eventurevent + 1)/2;
                else
                    trial_number = (EEG.epoch(1,i).eventurevent{1} + 1)/2;
                end
            end
        end
      
        if iscell(EEG.epoch(1,i).eventtype)
            event_code = EEG.epoch(1,i).eventtype{1};
        %Update: 1/22/2013
        %added to process event codes that are of string type
        elseif ischar(EEG.epoch(1,i).eventtype)
            event_code = str2double(EEG.epoch(1,i).eventtype);
        else
            event_code = EEG.epoch(1,i).eventtype(1);
        end
  
        if ischar(event_code)
            event_code = str2num(event_code);
        end
        
        if event_code == 11
            incong_idx(end+1,1)= i;
        end
        if event_code == 21
            cong_idx(end+1,1) = i;
        end
        if event_code ~=11 && event_code ~=21
            none_idx(end+1,1) = i;
        end
    end
    
    %Update: 1/22/2013
    %data integrity check
    if isempty(incong_idx)
        error('Failed to find first condition (incongruent) for subject %i probe',s)
    elseif isempty(cong_idx)
        error('Failed to find second condition (congruent) for subject %i probe',s)
    end
end

%% Update: 3/19/2013
%Block seperation
if boolUseBlockSeperation
    %loop through block list to find urevent range
    lstAcceptedRange = [];
    for iterator = 1:length(lstBlocks2Process)
        intBlockNum = lstBlocks2Process(iterator);
        idxStart = (intBlockNum - 1) * 82 + 1;
        idxEnd = idxStart + 82 - 1;
        lstAcceptedRange = [lstAcceptedRange,idxStart:idxEnd];
        clear intBlockNum idxStart idxEnd
    end
    lstAcceptedRange = sort(lstAcceptedRange);
    
    %get the urevents of all cong and incong events
    lstCongUrevents = [EEG.epoch(cong_idx).eventurevent];
    lstIncongUrevents = [EEG.epoch(incong_idx).eventurevent];
    
    %check if data type is cell
    if iscell(lstCongUrevents)
        lstCongUrevents = cell2mat(lstCongUrevents);
    end
    if iscell(lstIncongUrevents)
        lstIncongUrevents = cell2mat(lstIncongUrevents);
    end
    
    %if the urevents are bloated, use Glenn's algorithm to fix them
    %Note: subjects 20 and 110 have even urevents, using round to fix issue
    %of comparing integers to doubles
    if (max(lstCongUrevents) > 82 * 8 * 1.1)    %110% of 8 blocks of 82 trials, just to be safe
        lstCongUrevents = round((lstCongUrevents + 1)/2);
    end
    if (max(lstIncongUrevents) > 82 * 8 * 1.1)
        lstIncongUrevents = round((lstIncongUrevents + 1)/2);
    end
    
    %find the accepted events in the urevents lists
    lstAcceptedCongIndexes = find(ismember(lstCongUrevents,lstAcceptedRange));
    lstAcceptedIncongIndexes = find(ismember(lstIncongUrevents,lstAcceptedRange));
    
    %adjust the cong and incong indexes
    cong_idx = cong_idx(lstAcceptedCongIndexes);
    incong_idx = incong_idx(lstAcceptedIncongIndexes);
    
    %check if there is data left to process
    if isempty(incong_idx)
        error('Failed to find first condition (incongruent) for subject %i matching block range',s)
    elseif isempty(cong_idx)
        error('Failed to find second condition (congruent) for subject %i matching block range',s)
    end
end

%% start wavelet

%DEBUG: baseline subtraction exception
if DEBUG_DisableBaseline
    arrBaseline = NaN;
else
    arrBaseline = [-200 0];
end

%DEBUG: baseline subtraction method switching
if (DEBUG_FLAG_BaselineMethod == 0)
    strBaseline = 'off';
elseif (DEBUG_FLAG_BaselineMethod == 1)
    strBaseline = 'on';
elseif (DEBUG_FLAG_BaselineMethod == 2)
    strBaseline = 'full';
else
    error('Invalid baseline flag set')
end

%DEBUG: logrithmic scaling of frequency
if DEBUG_UseLogrithmicFrequencyScaling
    strFrequencyScaling = 'log';
else
    strFrequencyScaling = 'linear';
end

%-------------------------------------------------------------------------
  
  for n = 1:length(elecs)
        elec = elecs{n};

 if ischar(elec) == 1
        elec = str2num(elec);
 end

clc
disp(['Subject: ', num2str(s), ' | Electrode: ', num2str(n),'/',num2str(length(elecs))])

%DEBUG: warning messages
if DEBUG_DisableBaseline
    fprintf('Warning: baseline subtraction is disabled\n')
elseif (DEBUG_FLAG_BaselineMethod == 1)
    fprintf('Warning: baseline subtraction method is single-trial\n')
elseif (DEBUG_FLAG_BaselineMethod == 2)
    fprintf('Warning: baseline subtraction method is full\n')
end
if DEBUG_UseLogrithmicFrequencyScaling
    fprintf('Warning: frequency scaling is logrithmic\n')
end
if DEBUG_UsePermutation
    fprintf('Warning: permutation is on and will be saved\n')
end

        if DEBUG_UsePermutation
            [tmp_ersp itc powbase timesout freqs erspboot itcboot]=newtimef( {EEG.data(n,:,incong_idx) EEG.data(n,:,cong_idx)}, ...
            EEG.pnts,[EEG.xmin EEG.xmax]*1000, EEG.srate,c_level,'timesout',[times_define],'freqs',freqs, ...
            'elocs', EEG.chanlocs, 'baseline',arrBaseline,'trialbase',strBaseline,'freqscale',strFrequencyScaling,'alpha',DEBUG_intAlpha, ...
            'chaninfo', EEG.chaninfo,'padratio', 16, 'plotphase','off','plotersp','off','plotitc','off','verbose','on');
        
            %note: this is not fully implemented, the results are not saved
        else
            [tmp_ersp itc powbase timesout freqs]=newtimef( {EEG.data(n,:,incong_idx) EEG.data(n,:,cong_idx)}, ...
            EEG.pnts,[EEG.xmin EEG.xmax]*1000, EEG.srate,c_level,'timesout',[times_define],'freqs',freqs, ...
            'elocs', EEG.chanlocs, 'baseline',arrBaseline,'trialbase',strBaseline,'freqscale',strFrequencyScaling, ...
            'chaninfo', EEG.chaninfo,'padratio', 16, 'plotphase','off','plotersp','off','plotitc','off','verbose','on');
        end

        incong_ersp = tmp_ersp{1,1};
        cong_ersp   = tmp_ersp{1,2};
        sub_ersp    = tmp_ersp{1,3};
        
        if isequal(dose,'dose-1')
            dose = 'dose1'
        end


        eval_lead_incong     = ['ERSP_',group,'.',num2str(dose),'{',num2str(s),'}.incong.lead{',num2str(elec),'}=incong_ersp;']
        eval(eval_lead_incong);

        eval_lead_cong     = ['ERSP_',group,'.',num2str(dose),'{',num2str(s),'}.cong.lead{',num2str(elec),'}=cong_ersp;']
        eval(eval_lead_cong);

        eval_lead_contrast     = ['ERSP_',group,'.',num2str(dose),'{',num2str(s),'}.contrast.lead{',num2str(elec),'}=sub_ersp;']
        eval(eval_lead_contrast);
        
        if isequal(group,'ms') 
        if isequal(dose,'dose1')
            dose = 'dose-1'
        end
        end

        incong_itc = itc{1,1};
        cong_itc  = itc{1,2};
        sub_itc   = itc{1,3};
        
        if isequal(dose,'dose-1')
            dose = 'dose1'
        end


        eval_lead_incong     = ['ITC_',group,'.',num2str(dose),'{',num2str(s),'}.incong.lead{',num2str(elec),'}=incong_itc;']
        eval(eval_lead_incong);

        eval_lead_cong     = ['ITC_',group,'.',num2str(dose),'{',num2str(s),'}.cong.lead{',num2str(elec),'}=cong_itc;']
        eval(eval_lead_cong);

        eval_lead_contrast     = ['ITC_',group,'.',num2str(dose),'{',num2str(s),'}.contrast.lead{',num2str(elec),'}=sub_itc;']
        eval(eval_lead_contrast);
        
        if isequal(group,'ms') 
        if isequal(dose,'dose1')
            dose = 'dose-1'
        end
        end

  end

    
    cd(dirOutput)
    eval_save = ['save ERSP_',num2str(ends),' timesout freqs ERSP* History'];
    eval_save_itc = ['save ITC_',num2str(ends),' timesout freqs ITC* History'];
    eval(eval_save)
    eval(eval_save_itc)
end

if ~isempty(lstSubjectsNotFound)
    error('Failed to process these subjects: %s\n(all other subjects have finished processing)', num2str(lstSubjectsNotFound))
end