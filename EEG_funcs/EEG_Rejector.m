function [EEG, bad_CompOrChan]=EEG_Rejector(EEG_in,compType,params,varargin)
%EEG=EEG_Rejector(EEG,compType)
%Allowing viewing and rejecting bad channels or independent components
%Must load the respective file first
%file_dir: directory of the .set file during preprocessing
%compType: type of components
%1. rejecting bad channels
%2. rejecting bad independent components
%bad_compChan: channels or components rejected
%1. when rejecting bad channels, it is a cell array of channel names
%2. when rejecting bad independent components, it is a vector
%numComp: number of components to dispaly, required only when doing ICA
%flag='0';%default flag for prompting user's response
%no_plot: [0|1], turn off/on displaying plot, default turning on displaying
%plot (default no_plot=0)

%display prompts stored in a cell array
disp_lines={'channel(s)',...
    'Input as a cell array of strings {''A1'',...}:\n';...
    'components(s)',...
    'Input as a vector [1,2,3...]:\n'};

% flag: specify number of components to plot, only needed when viewing ICs
IND_numComp=find(strcmpi({'numComp'},varargin),1);
if ~isempty(IND_numComp)
    numComp=varargin{IND_numComp+1};
else
    numComp=15;%default numComp to be 15;
end
%flag: if needs to remove empty field
IND_PlotOff=find(strcmpi({'PlotOff'},varargin),1);
if ~isempty(IND_PlotOff)
    no_plot=varargin{IND_PlotOff+1};
else
    no_plot=0;%default: do not remove empty field
end



satisf=0;%a flag for whether or not user is satisfied with the result
counter=1;%count number of times one run is being done
while satisf==0
    %check if turning off plot
    if ~no_plot && numel(findobj(0,'type','figure')) <10
        switch compType%rejecting Channels
            case {1}%channels
                pop_eegplot(EEG_in,1,0,0)%display all channels for viewing
            case {2}%independent components
                %display components
                [~]=evalc('pop_prop(EEG_in,0,1:numComp,[])');
        end
    end
    %prompting users for behaviors
    disp(['Which ',disp_lines{compType,1},' would you like to reject?']);
    bad_CompOrChan=input(disp_lines{compType,2});
    if isempty(bad_CompOrChan)%exit if empty input, skip rejecting step
        EEG=EEG_in;%passing to output without modification
        return;
    end
    disp(['You are about to remove the following ',...
        disp_lines{compType,1},':']);
    disp(bad_CompOrChan);
    flag=input('Continue? (Y/N):','s');
    if strcmpi(flag,'y')
        switch compType
            case {1}%channels
                EEG = pop_select(EEG_in,'nochannel',bad_CompOrChan);
            case {2}%independent components
                EEG = pop_subcomp(EEG_in,bad_CompOrChan,0);
        end
    else
        continue;%start over
    end
    
    %dispaly plots again for double checking
    switch compType%rejecting Channels
        case {1}%channels
            pop_eegplot(EEG,1,0,0);%display all channels for viewing
        case {2}%independent components
            eegplot(EEG_in.data,'data2',EEG.data,'color',{'k'},...
                'title',['black = channel before rejection;'...
                ' red = after rejection -- eegplot()']);
            [~]=evalc(['[~,IND]=EEG_artifactRejection(EEG,2,',...
                '''thresh'',params.artThresh,''rejEpoch'',0)']);
            disp([num2str(length(IND)),' out of ',...
                num2str(length(EEG.epoch)), ' epoches will be rejected']);
    end
    
    %prompt user if satisfied with the result
    clear inquiry;
    disp('Satisfied with the results?');
    inquiry=input('A.Yes\nB.No,start over\n','s');
    switch upper(inquiry)
        case {'A','Y','YES'}
            satisf=1;%stop the loop and return the value
            close all;%close all the figures
        otherwise
            if compType==2 && counter>2
                sprintf('Increase the number of components displayed?\n');
                prompt_ans=input(['Enter the number of components to ',...
                    'display, or enter 0 to keep the same number of ',...
                    'components:\n']);
                if prompt_ans>0
                    numComp=prompt_ans;
                end
            end
            continue;%start over
    end
end
end