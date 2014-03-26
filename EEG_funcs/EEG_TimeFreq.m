function [ERSP,ITC,powbase,times_out,freqs_out]=...
    EEG_TimeFreq(DataSet,phase_in,Condition,subject,freqs_in,cycles,...
    varargin)
% DataSet is the RawDataSet of current phase
% Can sepcify if contrasting two conditions or not
% if turn off contrasting, ERSP and ITC returned will be a structure array
% that corresponds to each condition (should be in the order of the
% Condition table specified.
% phase_in: informatino of current phase
% Optional Inputs:
%       'TimeRes': time resolution, default 650
%       'contrast': turn on/off contrast, default 'on'


%inspect if Optional Inputs are altered
flag=InspectVarargin(varargin,{'contrast','TimeStep'},{1,5});

%get data
data_cell={DataSet.data.(Condition{1,1}),...
    DataSet.data.(Condition{2,1})};
eegparams=rmfield(DataSet,{'data','idx'});

%calculate times_in
times_in=DataSet.xmin*1000:flag.TimeStep:DataSet.xmax*1000;

switch flag.contrast
    case {1,'1','on','true','yes'}
        %get timeshift info
        timeshift=get_timeshift(phase_in,Condition(1:2,:));
            [ERSP,ITC,powbase,times_out,freqs_out]=...
                DO_newtimef(data_cell,eegparams,subject,phase_in,...
                times_in,freqs_in,cycles,timeshift);
    case {0,'0','off','false','no'}
        %place holding
        ERSP=struct();
        ITC=struct();
        powbase=cell(1,length(data_cell));
        times_out=cell(1,length(data_cell));
        freqs_out=cell(1,length(data_cell));
        %for each data set / condition
        for m = 1:length(data_cell)
            timeshift=get_timeshift(phase_in,Condition(m,:));
            data_mat=data_cell{m};
            [ERSP(m),ITC(m),powbase{m},times_out{m},freqs_out{m}]=...
                DO_newtimef(data_mat,eegparams,subject,phase_in,...
                times_in,freqs_in,cycles,timeshift);
        end
end
end

%% internal functions
function [ERSP,ITC,powbase,times_out,freqs_out]=DO_newtimef(...
    data_cell,eegparams,subject,phase_in,times_in,freqs_in,cycles,timeshift)
elecs = {eegparams.chanlocs.type};
%determine format of data
data_format=iscell(data_cell);
%do time frequency decomposition for each elecs
%display progress
h = waitbar(0,['Subject: ',subject,' | Phase: ',phase_in{1}]);
for n = 1:length(elecs)
    elec = elecs{n};
    if ischar(elec) == 1
        elec = str2num(elec);
    end
    %display progress
    waitbar(n/length(elecs),h,['Subject: ',subject,...
        ' | Phase: ',phase_in{1},' | Electrode: ', num2str(n)]);
    
    %reformat data_cell
    switch data_format
        case {1}%there are multiple data sets
            data_input=cellfun(@(x) x(n,:,:),data_cell,'Un',false);
        case {0}%only one matrix is specified
            data_input=data_cell(n,:,:);
    end 
    %may use one data set or two dataset
    [temp_ersp temp_itc powbase times_out freqs_out]=newtimef(...
        data_input, eegparams.pnts,[eegparams.xmin eegparams.xmax]*1000,...
        eegparams.srate, cycles,'timesout', times_in,'freqs',freqs_in,...
        'elocs', eegparams.chanlocs, 'baseline',timeshift*1000+[-200 0],...
        'chaninfo', eegparams.chaninfo,'plotphase','off',...
        'title',{[],[]},'plotersp','off','plotitc','off','verbose','off');
    
    %pass data to the ERSP and ITC data structure
    ERSP.incong.lead{elec}=temp_ersp{1,1};
    ERSP.cong.lead{elec}=temp_ersp{1,2};
    ERSP.contrast.lead{elec}=temp_ersp{1,3};
    
    ITC.incong.lead{elec}=temp_itc{1,1};
    ITC.cong.lead{elec}=temp_itc{1,2};
    ITC.contrast.lead{elec}=temp_itc{1,3};
end
delete(h);%remove progress bar
end


function timeshift=get_timeshift(phase_in,condition_in)
%get timeshift information from Condition; otherwise, no timeshift
if strcmpi(phase_in{5},'ON')
    timeshift=0;%if marker inserted at onset, no timeshift
else%if not, try to get timeshift
    try
        switch size(condition_in,1)
            case {1}
                timeshift=condition_in{3};
            case {2}
                timeshift=absExtreme(cell2mat(condition_in(:,3)),'max');
        end
        if isempty(timeshift)
            warning('No Time Shift Specified, use default timeshift 0');
            timeshift=0;%default: no timeshift
        end
    catch
        warning('No Time Shift Specified, use default timeshift 0');
        timeshift=0;%default: no timeshift
    end
end

end