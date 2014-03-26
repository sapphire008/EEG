%Headplot pipeline
%Author: Julian Y. Cheng
%2/21/2013
%
%This script is an adaption of srPlotHeadplots and srPermutation_headplot
%to combine the functionality into a single processing pipeline as well as
%make use of the new data structure.
%
%Note: the permutation section of this script makes use of the new
%permutation_test function by Scott; it requires that said section be run
%on a machine equiped with a NVIDIA GPU for CUDA
%Note: only ERSP's are processed, ITC's are ignored. In the future if ITC
%is needed, used the source scripts that each code cell is adapted from as
%reference to reconstruct the ITC processing code.

%% Single plot
%
%This code cell processes a single input file per user selection. Used for
%debugging purposes. Adapted from srPlotSpectrograms
%
%Changelog:
%   3/28/2013:  Added processing of History variable if it exists. If
%               multiple files exist for the same data type (i.e. if
%               subject groups were processed in different batches), then
%               the History from the first file will be used.

clc,clear

%user-defined vars

%general
arrFrequencyRange = [4 8];  %Hz
arrTimeInterval = [0 250];  %ms
arrRotation = [0 90];       %3D rotation for 3D topographic plots;
                            %[0     90]: top down
                            %[-100  0 ]: left lateral
                            %[100   0 ]: right lateral
intRejectLimit = 17;        %maximum number of subjects that can have empty data per electrode

%subject exclusion
boolUseSubjectExclusion = false;
lstExcludedSubjects = [103,108,110,111,115];   %will not process these subjects, use to limit the MC/PC dataset to MC or PC only

%task-specific vars

%general
dirRoot = '/nfs/erp-modaf/mc/MC-epoched-Feb11-Glenn/';
strInputFolder = 'Wavelet results/';                                    %combined with dirRoot and phase string to produce path to input folders
strSubfolder = 'Processed matrices/';                                   %the folder to search for the input MAT files
strSplineFile = 'Headplots/HDPLT_SPLINE.spl';                           %spline file will be saved here, combined with dirRoot and strPhase
pathCleanEEGSet = '/nfs/erp-modaf/mc/ERP/Generic_set_1trial.set';       %needs to have complete chanloc information for all channels
pathElectrodeMap = '/nfs/erp-modaf/elec_files/wholehead_elecs.mat';     %defines the regions
intElectrodes = 128;                                                    %total electrode count, used to construct electrode map of all brain
lstMode = {'2D topographic','3D topographic'};
lstPhase = {'Cue','Probe'};

%comparisons
%note: this specifies which comparisons are graphable, by defining the
%required files here. Note that for Qmoney conditions, the order of the
%files define the order of subtraction; i.e. file1 - file2
structComparisons(1).name = 'Drug green';
structComparisons(1).files = {'Drug_cong'}; 
structComparisons(2).name = 'Drug red';
structComparisons(2).files = {'Drug_incong'}; 
structComparisons(3).name = 'Drug contrast';
structComparisons(3).files = {'Drug_contrast'}; 
structComparisons(4).name = 'Plac green';
structComparisons(4).files = {'Plac_cong'}; 
structComparisons(5).name = 'Plac red';
structComparisons(5).files = {'Plac_incong'}; 
structComparisons(6).name = 'Plac contrast';
structComparisons(6).files = {'Plac_contrast'}; 
structComparisons(7).name = 'Drug green - Plac green';
structComparisons(7).files = {'Drug_cong','Plac_cong'}; 
structComparisons(8).name = 'Drug red - Plac red';
structComparisons(8).files = {'Drug_incong','Plac_incong'}; 
structComparisons(9).name = 'Qmoney';
structComparisons(9).files = {'Drug_contrast','Plac_contrast'}; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%begin algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%regions definitions-------------------------------------------------------

structRegions(1).name = 'all-Head';
structRegions(1).electrodes = 1:intElectrodes;

load(pathElectrodeMap)

structRegions(2).name = 'all-PFC';
structRegions(2).electrodes = cat(2, wholehead_elecs{1,2},wholehead_elecs{2,2}, wholehead_elecs{3,2});
structRegions(3).name = 'left-PFC';
structRegions(3).electrodes = wholehead_elecs{2,2};
structRegions(4).name = 'mid-PFC';
structRegions(4).electrodes = wholehead_elecs{1,2};
structRegions(5).name = 'right-PFC';
structRegions(5).electrodes = wholehead_elecs{3,2};
structRegions(6).name = 'all-Central';
structRegions(6).electrodes = cat(2, wholehead_elecs{8,2},wholehead_elecs{4,2}, wholehead_elecs{9,2});
structRegions(7).name = 'left-Central';
structRegions(7).electrodes = wholehead_elecs{8,2};
structRegions(8).name = 'mid-Central';
structRegions(8).electrodes = wholehead_elecs{4,2};
structRegions(9).name = 'right-Central';
structRegions(9).electrodes = wholehead_elecs{9,2};
structRegions(10).name = 'all-Post';
structRegions(10).electrodes = cat(2, wholehead_elecs{6,2},wholehead_elecs{5,2}, wholehead_elecs{7,2});
structRegions(11).name = 'left-Post';
structRegions(11).electrodes = wholehead_elecs{6,2};
structRegions(12).name = 'mid-Post';
structRegions(12).electrodes = wholehead_elecs{5,2};
structRegions(13).name = 'right-Post';
structRegions(13).electrodes = wholehead_elecs{7,2};

clear wholehead_elecs

%begin ui algorithm--------------------------------------------------------

%prompt for mode
[idxSelection boolOK] = listdlg('Name','Mode','PromptString','Select operation mode:', ...
                                'ListString',lstMode,'ListSize',[200 300],'SelectionMode','single');
if boolOK
    flagMode = idxSelection;
else
    error('Operation aborted')
end
clear idxSelection boolOK

%prompt for phase

%show dialog
[idxSelection boolOK] = listdlg('Name','Phase','PromptString','Select phase:', ...
                                'ListString',lstPhase,'ListSize',[200 300],'SelectionMode','single');
if boolOK
    strPhase = lstPhase{idxSelection};
else
    error('Operation aborted')
end
clear lstPhase idxSelection boolOK

%prompt for analysis folder

%generate path to input folders
dirInputFolder = [dirRoot,strPhase,'/',strInputFolder];
pathSpline = [dirRoot,strPhase,'/',strSplineFile];

%find all analysis folders
lstInfolders = dir(dirInputFolder);
lstInfolders = {lstInfolders(3:end).name};  %the indexing is to get rid of "." and ".."

%show dialog
[idxSelection boolOK] = listdlg('Name','Analysis','PromptString','Select analysis folder:', ...
                                'ListString',lstInfolders,'ListSize',[300 300],'SelectionMode','single');
if boolOK
    strInfolder = lstInfolders{idxSelection};
else
    error('Operation aborted')
end
clear lstInfolders idxSelection boolOK

%sanity check: make sure subfolder exists
dirAnalysis = [dirInputFolder,strInfolder,'/',strSubfolder];
if ~exist(dirAnalysis,'dir')
    error('Failed to find path to processed matrices: %s',dirAnalysis)
end

%prompt for data type

%find all processed matrices
lstInfiles = dir(dirAnalysis);          %this variable is not destroyed in ui algorithm
lstInfiles = {lstInfiles(3:end).name};
idxDelete = find(~cellfun(@isempty,strfind(lstInfiles,'ITC')));
lstInfiles(idxDelete) = [];
clear idxDeletew

%prompt for comparisons to graph

%show dialog
[idxSelection boolOK] = listdlg('Name','Comparison','PromptString','Select comparison(s):', ...
                                'ListString',{structComparisons.name},'ListSize',[300 300],'SelectionMode','multiple');
if boolOK
    structComparisons = structComparisons(idxSelection);
else
    error('Operation aborted')
end
clear idxSelection boolOK

%prompt for regions to graph

%show dialog
[idxSelection boolOK] = listdlg('Name','Region','PromptString','Select region(s):', ...
                                'ListString',{structRegions.name},'ListSize',[200 300],'SelectionMode','multiple');
if boolOK
    structRegions = structRegions(idxSelection);
else
    error('Operation aborted')
end
clear idxSelection boolOK

%begin processing algorithm------------------------------------------------

% %get chanlocs from a complete .set file
[~] = evalc('structEEG = pop_loadset(pathCleanEEGSet);');   %silence the function output
structChanlocs = structEEG.chanlocs;

%loop through all comparisons
for myComparison = structComparisons
    
    %loop through all required files
    lstFileIndexes = cell(1,length(myComparison.files));
    for i = 1:length(myComparison.files)
        strFile = myComparison.files{i};
        
        %save indexes of required file found among infiles
        lstFileIndexes{i} = find(~cellfun(@isempty,strfind(lstInfiles,strFile)));
    end
    clear strFile i
    
    %make sure all required files have the same count (at least 1)
    arrFileCount = cellfun(@length,lstFileIndexes);
    if ~all(arrFileCount) || ~all(arrFileCount(:) == arrFileCount(1))
        %if not all numbers are non-zero OR not all numbers are the same
        fprintf('Warning: required file matching failed in %s, skipped\n',myComparison.name)
        continue
    end
    clear arrFileCount
    
    %load all files
    %note: this is coded for maximum flexibility, in practice only 1 file
    %will be loaded for Drug/Plac comparisons and 2 files for Qmoney
    %comparisons (if both ERSP and ITC are graphed, the numbers will double)
    lstFileIndexes = cell2mat(lstFileIndexes(:)');  %discard sorting by search string, instead just load all files like in pipePermutationTest
    for i = 1:length(lstFileIndexes)
        load([dirAnalysis,lstInfiles{lstFileIndexes(i)}])
        
        %save into data structure
        dataWavelet(i) = Data;
        if exist('History','var')
            dataHistory(i) = History;
        end
        
        clear Data History
    end
    clear lstFileIndexes
    
    %loop through all regions
    for myRegion = structRegions
        %get the electrode list
        lstElectrodes = myRegion.electrodes;
        
        %process both ERSP and ITC
        %note: this code is not changed because if ITC doesn't exist it is
        %automatically skipped
        for cellDatatype = {'ERSP','ITC'}
            strDatatype = cell2mat(cellDatatype);
            
            %find indexes of current data type
            idxDatatype = find(cellfun(@(x) strcmp(x,strDatatype),{dataWavelet.type}));
            
            %check if data is present
            if isempty(idxDatatype);
                %no data is present for current data type, skip
                continue
            end
            
            %check if there is more than one set of data to process
            if (length(idxDatatype) > 1)
                %parse dosage and condition information from myComparison.files
                strSplit = regexp(myComparison.files{1},'_','split');
                strDosage1 = strSplit{1};
                strCondition1 = strSplit{2};
                strSplit = regexp(myComparison.files{2},'_','split');
                strDosage2 = strSplit{1};
                strCondition2 = strSplit{2};
                
                %find where the conditions are stored in the dataWavelet structure
                idxDosage = find(cellfun(@(x) strcmp(x,strDosage1),{dataWavelet.field}));
                idxCondition = find(cellfun(@(x) strcmp(x,strCondition1),{dataWavelet.subfield}));
                idxGroup1 = intersect(idxDosage,idxCondition);
                idxGroup1 = intersect(idxDatatype,idxGroup1);
                idxDosage = find(cellfun(@(x) strcmp(x,strDosage2),{dataWavelet.field}));
                idxCondition = find(cellfun(@(x) strcmp(x,strCondition2),{dataWavelet.subfield}));
                idxGroup2 = intersect(idxDosage,idxCondition);
                idxGroup2 = intersect(idxDatatype,idxGroup2);
                
                %perform subtraction, and modify structure of Data
                myData = dataWavelet(idxGroup1);
                myData.field = myComparison.name;
                myData.subfield = '';   %this is ignored
                myData.matrix = dataWavelet(idxGroup1).matrix - dataWavelet(idxGroup2).matrix;
                            
                %build source string
                if exist('dataHistory','var')
                    strSource1_name = dataHistory(idxGroup1).source_name;
                    strSource1_id = dataHistory(idxGroup1).source_id;
                    strSource2_name = dataHistory(idxGroup2).source_name;
                    strSource2_id = dataHistory(idxGroup2).source_id;
                    if (strcmp(strSource1_name,strSource2_name) && strcmp(strSource1_id,strSource2_id))
                        strSource = ['Source: ',strrep(strSource1_name,'_','\_'),' (',strSource1_id,')'];
                    else
                        strSource = {['Source: ',strrep(strSource1_name,'_','\_'),' (',strSource1_id,')'],[strrep(strSource2_name,'_','\_'),' (',strSource2_id,')']};
                    end
                else
                    strSource = ['Source: ',strrep(myData.source,'_','\_')];
                end
                
                clear strDosage1 strDosage2 strCondition1 strCondition2
                clear idxDosage idxCondition idxGroup1 idxGroup2
                clear strSource1_name strSource1_id strSource2_name strSource2_id
            else
                myData = dataWavelet(idxDatatype);
                %build source string
                if exist('dataHistory','var')
                    strSource_name = dataHistory(idxDatatype).source_name;
                    strSource_id = dataHistory(idxDatatype).source_id;
                    strSource = ['Source: ',strrep(strSource_name,'_','\_'),' (',strSource_id,')'];
                else
                    strSource = ['Source: ',strrep(myData.source,'_','\_')];
                end
                clear strSource_name strSource_id
            end
                
            %check if subject exclusion is enabled
            if boolUseSubjectExclusion
                idxDelete = find(ismember(myData.subject_list,lstExcludedSubjects));
                myData.matrix(idxDelete,:,:,:) = [];
                
                clear idxDelete
            end
            
            %extract data by frequency
            arrFrequencyIndexes = find(myData.frequency >= arrFrequencyRange(1),1):find(myData.frequency <= arrFrequencyRange(2),1,'last');
            myData.matrix = myData.matrix(:,:,arrFrequencyIndexes,:);
            clear arrFrequencyIndexes
            
            %extract data by time
            arrTimesIndexes = find(myData.time >= arrTimeInterval(1),1):find(myData.time <= arrTimeInterval(2),1,'last');
            myData.matrix = myData.matrix(:,:,:,arrTimesIndexes);
            clear arrTimesIndexes
            
            %average down the frequency <subject X electrode X time>
            myData.matrix = squeeze(nanmean(myData.matrix,3));
            
            %average down the time <subject X electrode>
            myData.matrix = squeeze(nanmean(myData.matrix,3));
            
            %find electrodes to exclude
            idxElectrodes2Remove = [];
            for i = 1:size(myData.matrix,2)
                intExcludedSubjects = sum(isnan(myData.matrix(:,i)));
                if (intExcludedSubjects > intRejectLimit)
                    idxElectrodes2Remove(end+1) = i;
                    fprintf('Removed electrode %i for having %i excluded subjects\n',i,intExcludedSubjects);
                end
            end
            idxElectrodes2Remove = union(idxElectrodes2Remove,setdiff(1:intElectrodes,lstElectrodes));
            lstElectrodes2Graph = setdiff(1:intElectrodes,idxElectrodes2Remove);
            
            %average down the subject dimension
            myData.matrix = nanmean(myData.matrix,1);
            
            %graph the result
            switch flagMode
                case 1  %2D
                    figure,colorbar
                    h = topoplot(myData.matrix,structChanlocs,'electrodes','on','emarker',{'.','k',5,1},'plotchans',lstElectrodes2Graph);
                case 2  %3D
                    evalc('structEEG = pop_select(structEEG,''nochannel'',idxElectrodes2Remove)');
                    evalc('headplot(''setup'',structEEG.chanlocs,pathSpline)'); %create spline file
                    h = headplot(myData.matrix(lstElectrodes2Graph),pathSpline,'view',arrRotation,'cbar',0,'electrodes','on','verbose','off');
                otherwise
                    error('Invalid graphing mode selection')
            end
            lblText = axes('Units','Normal','Position',[0.35 -0.86 .85 .85],'Visible','off');
            set(get(lblText,'Title'),'Visible','on'), title(strSource,'FontSize',8);
            
            clear myData idxDatatype strDatatype lblText idxColorScale lstElectrodes2Graph idxElectrodes2Remove strSource
        end
        clear lstElectrodes
    end
end

%% Multi plot
%
%This code cell is used to graph all plots from a single input folder. Uses
%UI based switching for phase, analysis folder selection, regions, and 
%output folder specification. Adapted from srPlotSpectrograms
%
%Changelog:
%   3/28/2013:  Added processing of History variable if it exists. If
%               multiple files exist for the same data type (i.e. if
%               subject groups were processed in different batches), then
%               the History from the first file will be used. Also, 3D
%               rotation graphing is adopted back from SC.

clc,close,clear

%user-defined vars

%general
arrFrequencyRange = [15 30];  %Hz
arrTimeInterval = [0 200];  %ms
intRejectLimit = 17;        %maximum number of subjects that can have empty data per electrode
boolPublicationMode = false; %this disables titles and source labeling in the graphs

%subject exclusion
boolUseSubjectExclusion = false;
lstExcludedSubjects = [103,108,110,111,115];   %will not process these subjects, use to limit the MC/PC dataset to MC or PC only

%task-specific vars

%general
dirRoot = '/nfs/erp-modaf/mc/MC-epoched-Feb11-Glenn/';
strInputFolder = 'Wavelet results/';                                    %combined with dirRoot and phase string to produce path to input folders
strSubfolder = 'Processed matrices/';                                   %the folder to search for the input MAT files
strOutputFolder = 'Headplots/Figures/';
strPublicationFolder = 'Publication/';                                  %this is the subfolder that would be created within the normal output directory in publication mode
strSplineFile = 'Headplots/HDPLT_SPLINE.spl';                           %spline file will be saved here, combined with dirRoot and strPhase
strFilePrefix = 'HDPLT';
pathCleanEEGSet = '/nfs/erp-modaf/mc/ERP/Generic_set_1trial.set';       %needs to have complete chanloc information for all channels
intElectrodes = 128;                                                    %total electrode count, used to construct electrode map of all brain
lstMode = {'2D topographic','3D topographic'};
lstPhase = {'Cue','Probe'};

%comparisons
%note: this specifies which comparisons are graphable, by defining the
%required files here. Note that for Qmoney conditions, the order of the
%files define the order of subtraction; i.e. file1 - file2
structComparisons(1).name = 'Drug green';
structComparisons(1).files = {'Drug_cong'}; 
structComparisons(2).name = 'Drug red';
structComparisons(2).files = {'Drug_incong'}; 
structComparisons(3).name = 'Drug contrast';
structComparisons(3).files = {'Drug_contrast'}; 
structComparisons(4).name = 'Plac green';
structComparisons(4).files = {'Plac_cong'}; 
structComparisons(5).name = 'Plac red';
structComparisons(5).files = {'Plac_incong'}; 
structComparisons(6).name = 'Plac contrast';
structComparisons(6).files = {'Plac_contrast'}; 
structComparisons(7).name = 'Drug green - Plac green';
structComparisons(7).files = {'Drug_cong','Plac_cong'}; 
structComparisons(8).name = 'Drug red - Plac red';
structComparisons(8).files = {'Drug_incong','Plac_incong'}; 
structComparisons(9).name = 'Qmoney';
structComparisons(9).files = {'Drug_contrast','Plac_contrast'}; 

%3D rotation views
%note: ignored in 2D mode
structRotations(1).name = 'top down';
structRotations(1).array = [0,90];
structRotations(1).postfix = '';
structRotations(2).name = 'left lateral';
structRotations(2).array = [-100,0];
structRotations(2).postfix = '_left';
structRotations(3).name = 'right lateral';
structRotations(3).array = [100,0];
structRotations(3).postfix = '_right';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%begin algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%regions definitions-------------------------------------------------------

structRegions(1).name = 'all-Head';
structRegions(1).electrodes = 1:intElectrodes;

%begin ui algorithm--------------------------------------------------------

%prompt for mode
[idxSelection boolOK] = listdlg('Name','Mode','PromptString','Select operation mode:', ...
                                'ListString',lstMode,'ListSize',[200 300],'SelectionMode','single');
if boolOK
    flagMode = idxSelection;
else
    error('Operation aborted')
end
clear idxSelection boolOK

%prompt for phase

%show dialog
[idxSelection boolOK] = listdlg('Name','Phase','PromptString','Select phase:', ...
                                'ListString',lstPhase,'ListSize',[200 300],'SelectionMode','single');
if boolOK
    strPhase = lstPhase{idxSelection};
else
    error('Operation aborted')
end
clear lstPhase idxSelection boolOK

%prompt for analysis folder

%generate path to input folders
dirInputFolder = [dirRoot,strPhase,'/',strInputFolder];
pathSpline = [dirRoot,strPhase,'/',strSplineFile];

%find all analysis folders
lstInfolders = dir(dirInputFolder);
lstInfolders = {lstInfolders(3:end).name};  %the indexing is to get rid of "." and ".."

%show dialog
[idxSelection boolOK] = listdlg('Name','Analysis','PromptString','Select analysis folder:', ...
                                'ListString',lstInfolders,'ListSize',[300 300],'SelectionMode','single');
if boolOK
    strInfolder = lstInfolders{idxSelection};
else
    error('Operation aborted')
end
clear lstInfolders idxSelection boolOK

%sanity check: make sure subfolder exists
dirAnalysis = [dirInputFolder,strInfolder,'/',strSubfolder];
if ~exist(dirAnalysis,'dir')
    error('Failed to find path to processed matrices: %s',dirAnalysis)
end

%find all processed matrices
lstInfiles = dir(dirAnalysis);          %this variable is not destroyed in ui algorithm
lstInfiles = {lstInfiles(3:end).name};
idxDelete = find(~cellfun(@isempty,strfind(lstInfiles,'ITC')));
lstInfiles(idxDelete) = [];
clear idxDelete

%prompt for output folder

%create selection list
lstOutput = {'{Frequency min/max} {Time min/max} format','Custom'};

%show dialog
[idxSelection boolOK] = listdlg('Name','Output directory','PromptString','Select output directory name setting:', ...
                                'ListString',lstOutput,'ListSize',[200 300],'SelectionMode','single');
if boolOK
    switch idxSelection
        case 1  %frequency-time format
            strBuilder = [num2str(arrFrequencyRange(1)),'-',num2str(arrFrequencyRange(2)),'Hz '];
            strBuilder = [strBuilder,num2str(arrTimeInterval(1)),'-',num2str(arrTimeInterval(2)),'ms'];
            strOutfolder = strBuilder;
            clear strBuilder
        case 2  %custom
            %setup input textbox
            paramPrompt = {'Enter custom folder name:'};
            paramName = 'Input custom folder name';
            paramLines = 1;
            paramDefault = {''};
            paramOptions.Resize = 'on';
            paramOptions.WindowStyle = 'modal';
            paramOptions.Interpreter = 'none';
            
            %show dialog
            cellResponse = inputdlg(paramPrompt,paramName,paramLines,paramDefault,paramOptions);
            
            if isempty(cellResponse)
                error('Operation aborted')
            elseif isempty(cell2mat(cellResponse))
                error('Input dialog failure: folder name cannot be empty')
            else
                strOutfolder = strrep(cell2mat(cellResponse),'/','');   %strrep to remove / characters; this is manually added later
            end
            clear paramPrompt paramName paramLines paramDefault paramOptions cellResponse
        otherwise
            error('Impossible selection reached for output directory setting')
    end
else
    error('Operation aborted')
end
clear idxSelection boolOK lstOutput

%prompt for overwrite protection

if exist([dirRoot,strPhase,'/',strOutputFolder,strOutfolder,'/'],'dir')
    %show dialog
    strResponse = questdlg('Folder already exists, overwrite? (NO will append numbers as necessary to folder name)','Overwrite');
    
    switch strResponse
        case 'Yes'
            dirOutput = [dirRoot,strPhase,'/',strOutputFolder,strOutfolder,'/'];
        case 'No'
            intCount = 2;
            while exist([dirRoot,strPhase,'/',strOutputFolder,strOutfolder,'-',num2str(intCount),'/'],'dir')
                intCount = intCount +1;
            end
            dirOutput = [dirRoot,strPhase,'/',strOutputFolder,strOutfolder,'-',num2str(intCount),'/'];
            mkdir(dirOutput);
            clear intCount
        otherwise
            error('Operation aborted')
    end
    clear strResponse
else
    dirOutput = [dirRoot,strPhase,'/',strOutputFolder,strOutfolder,'/'];
    mkdir(dirOutput);
end

if boolPublicationMode
    dirOutput = [dirOutput,strPublicationFolder];
    if ~exist(dirOutput,'dir')
        mkdir(dirOutput);
    end
end

%begin processing algorithm------------------------------------------------

%load all files
h = waitbar(0,'Loading... ');
for i = 0:length(lstInfiles)
    if (i == 0)
        %special case used to load clean EEG set file
        waitbar(i/(length(lstInfiles)+1),h,'Loading... Clean EEG set');
        [~] = evalc('structEEG = pop_loadset(pathCleanEEGSet);');   %silence the function output
        structChanlocs = structEEG.chanlocs;
        waitbar((i+1)/(length(lstInfiles)+1),h);
        continue
    end
    
    waitbar(i/(length(lstInfiles)+1),h,['Loading... ',strrep(lstInfiles{i},'_','\_')]);
    
    load([dirAnalysis,lstInfiles{i}])

    %save into data structure
    dataWavelet(i) = Data;
    if exist('History','var')
        dataHistory(i) = History;
    end
    
    waitbar((i+1)/(length(lstInfiles)+1),h);

    clear Data History
end
delete(h);

%loop through all comparisons
figure
lstGraphedElectrodes = [];  %this is used to avoid recompiling 3D meshes if the electrodes graphed did not change due to exclusion criteria
for myComparison = structComparisons
    
    %loop through all regions
    for myRegion = structRegions
        %get the electrode list
        lstElectrodes = myRegion.electrodes;
        
        %process both ERSP and ITC
        %note: this code is not changed because if ITC doesn't exist it is
        %automatically skipped
        for cellDatatype = {'ERSP','ITC'}
            strDatatype = cell2mat(cellDatatype);
            
            %find indexes of current data type
            idxDatatype = find(cellfun(@(x) strcmp(x,strDatatype),{dataWavelet.type}));
            
            %check if data is present
            if isempty(idxDatatype);
                %no data is present for current data type, skip
                continue
            end
            
            %find the data set from the collection dataWavelet
            %check if manual subtraction is needed (for Qmoney comparisons)
            if (length(myComparison.files) > 1)
                %parse dosage and condition information from myComparison.files
                strSplit = regexp(myComparison.files{1},'_','split');
                strDosage1 = strSplit{1};
                strCondition1 = strSplit{2};
                strSplit = regexp(myComparison.files{2},'_','split');
                strDosage2 = strSplit{1};
                strCondition2 = strSplit{2};
                
                %find where the conditions are stored in the dataWavelet structure
                idxDosage = find(cellfun(@(x) strcmp(x,strDosage1),{dataWavelet.field}));
                idxCondition = find(cellfun(@(x) strcmp(x,strCondition1),{dataWavelet.subfield}));
                idxGroup1 = intersect(idxDosage,idxCondition);
                idxGroup1 = intersect(idxDatatype,idxGroup1);
                idxDosage = find(cellfun(@(x) strcmp(x,strDosage2),{dataWavelet.field}));
                idxCondition = find(cellfun(@(x) strcmp(x,strCondition2),{dataWavelet.subfield}));
                idxGroup2 = intersect(idxDosage,idxCondition);
                idxGroup2 = intersect(idxDatatype,idxGroup2);
                
                %perform subtraction, and modify structure of Data
                myData = dataWavelet(idxGroup1);
                myData.field = myComparison.name;
                myData.subfield = '';   %this is ignored
                myData.matrix = dataWavelet(idxGroup1).matrix - dataWavelet(idxGroup2).matrix;
                            
                %build source string
                if exist('dataHistory','var')
                    strSource1_name = dataHistory(idxGroup1).source_name;
                    strSource1_id = dataHistory(idxGroup1).source_id;
                    strSource2_name = dataHistory(idxGroup2).source_name;
                    strSource2_id = dataHistory(idxGroup2).source_id;
                    if (strcmp(strSource1_name,strSource2_name) && strcmp(strSource1_id,strSource2_id))
                        strSource = ['Source: ',strrep(strSource1_name,'_','\_'),' (',strSource1_id,')'];
                    else
                        strSource = {['Source: ',strrep(strSource1_name,'_','\_'),' (',strSource1_id,')'],[strrep(strSource2_name,'_','\_'),' (',strSource2_id,')']};
                    end
                else
                    strSource = ['Source: ',strrep(myData.source,'_','\_')];
                end
                
                clear strDosage1 strDosage2 strCondition1 strCondition2
                clear idxDosage idxCondition idxGroup1 idxGroup2
                clear strSource1_name strSource1_id strSource2_name strSource2_id
            else
                %parse dosage and condition information from myComparison.files
                strSplit = regexp(cell2mat(myComparison.files),'_','split');
                strDosage = strSplit{1};
                strCondition = strSplit{2};
                
                %find where the conditions are stored in the dataWavelet structure
                idxDosage = find(cellfun(@(x) strcmp(x,strDosage),{dataWavelet.field}));
                idxCondition = find(cellfun(@(x) strcmp(x,strCondition),{dataWavelet.subfield}));
                idxGroup = intersect(idxDosage,idxCondition);
                idxGroup = intersect(idxDatatype,idxGroup);
                
                %get the data
                myData = dataWavelet(idxGroup);
                
                %build source string
                if exist('dataHistory','var')
                    strSource_name = dataHistory(idxGroup).source_name;
                    strSource_id = dataHistory(idxGroup).source_id;
                    strSource = ['Source: ',strrep(strSource_name,'_','\_'),' (',strSource_id,')'];
                else
                    strSource = ['Source: ',strrep(myData.source,'_','\_')];
                end
                
                clear strDosage strCondition idxDosage idxCondition idxGroup
                clear strSource_name strSource_id
            end
                
            %check if subject exclusion is enabled
            if boolUseSubjectExclusion
                idxDelete = find(ismember(myData.subject_list,lstExcludedSubjects));
                myData.matrix(idxDelete,:,:,:) = [];
                
                clear idxDelete
            end
            
            %extract data by frequency
            arrFrequencyIndexes = find(myData.frequency >= arrFrequencyRange(1),1):find(myData.frequency <= arrFrequencyRange(2),1,'last');
            myData.matrix = myData.matrix(:,:,arrFrequencyIndexes,:);
            clear arrFrequencyIndexes
            
            %extract data by time
            arrTimesIndexes = find(myData.time >= arrTimeInterval(1),1):find(myData.time <= arrTimeInterval(2),1,'last');
            myData.matrix = myData.matrix(:,:,:,arrTimesIndexes);
            clear arrTimesIndexes
            
            %average down the frequency <subject X electrode X time>
            myData.matrix = squeeze(nanmean(myData.matrix,3));
            
            %average down the time <subject X electrode>
            myData.matrix = squeeze(nanmean(myData.matrix,3));
            
            %find electrodes to exclude
            idxElectrodes2Remove = [];
            for i = 1:size(myData.matrix,2)
                intExcludedSubjects = sum(isnan(myData.matrix(:,i)));
                if (intExcludedSubjects > intRejectLimit)
                    idxElectrodes2Remove(end+1) = i;
                    fprintf('Removed electrode %i for having %i excluded subjects in %s\n',i,intExcludedSubjects,myComparison.name);
                end
            end
            lstElectrodes2Graph = lstElectrodes;
            lstElectrodes2Graph(idxElectrodes2Remove) = [];
            
            %average down the subject dimension
            myData.matrix = nanmean(myData.matrix,1);
            
            %build filename
            strBuilder = [strFilePrefix,'_'];
            strBuilder = [strBuilder,strDatatype,'_'];
            strBuilder = [strBuilder,myComparison.name,'_'];
            strBuilder = [strBuilder,num2str(arrFrequencyRange(1)),'-',num2str(arrFrequencyRange(2)),'Hz_'];
            strBuilder = [strBuilder,num2str(arrTimeInterval(1)),'-',num2str(arrTimeInterval(2)),'ms_'];
            switch flagMode
                case 1  %2D
                    strBuilder = [strBuilder,'2D'];
                case 2  %3D
                    %moved this to rotation for loop below
                    %strBuilder = [strBuilder,'3D',strPostfix];
                otherwise
                    error('Invalid graphing mode selection')
            end
            
            %graph the result
            clf;
            switch flagMode
                case 1  %2D
                    colorbar;
                    h = topoplot(myData.matrix,structChanlocs,'electrodes','on','emarker',{'.','k',5,1},'plotchans',lstElectrodes2Graph);
                    if ~boolPublicationMode
                        lblText = axes('Units','Normal','Position',[0.35 -0.86 .85 .85],'Visible','off');
                        set(get(lblText,'Title'),'Visible','on'), title(strSource,'FontSize',6);
                    end
                    saveas(gcf,[dirOutput,strBuilder,'.tif'],'tiff');
                case 2  %3D
                    %only setup spline if the list of electrodes to graphchanged
                    if ~isequal(lstElectrodes2Graph,lstGraphedElectrodes)
                        evalc('headplot(''setup'',structChanlocs(lstElectrodes2Graph),pathSpline)'); %create spline file
                        lstGraphedElectrodes = lstElectrodes2Graph;
                    end
                    h = headplot(myData.matrix(lstElectrodes2Graph),pathSpline,'cbar',0,'electrodes','on','verbose','off');
                    if ~boolPublicationMode
                        lblText = axes('Units','Normal','Position',[0.35 -0.86 .85 .85],'Visible','off');
                        set(get(lblText,'Title'),'Visible','on'), title(strSource,'FontSize',6);
                        axes(h);
                    end
                    %process all rotation views
                    for i = 1:length(structRotations)
                        view(structRotations(i).array);
                        saveas(gcf,[dirOutput,strBuilder,'3D',structRotations(i).postfix,'.tif'],'tiff');
                    end
                otherwise
                    error('Invalid graphing mode selection')
            end
            
            %old algorithm, archived in case new code has bugs
%             %graph the result
%             clf;
%             switch flagMode
%                 case 1  %2D
%                     colorbar;
%                     h = topoplot(myData.matrix,structChanlocs,'electrodes','on','emarker',{'.','k',5,1},'plotchans',lstElectrodes2Graph);
%                 case 2  %3D
%                     evalc('structEEG = pop_select(structEEG,''nochannel'',idxElectrodes2Remove)');
%                     evalc('headplot(''setup'',structEEG.chanlocs,pathSpline)'); %create spline file
%                     h = headplot(myData.matrix(lstElectrodes2Graph),pathSpline,'view',arrRotation,'cbar',0,'electrodes','on','verbose','off');
%                 otherwise
%                     error('Invalid graphing mode selection')
%             end
%             if ~boolPublicationMode
%                 lblText = axes('Units','Normal','Position',[0.3 -0.86 .85 .85],'Visible','off');
%                 set(get(lblText,'Title'),'Visible','on'), title(strSource,'FontSize',6);
%             end
%             saveas(gcf,[dirOutput,strBuilder,'.tif'],'tiff');
            
            clear myData idxDatatype strDatatype lblText idxColorScale strBuilder lstElectrodes2Graph idxElectrodes2Remove
        end
        clear lstElectrodes
    end
end

close,beep

%% Permutation
%
%This code cell is an adaption of pipePermutationTest
%
%Note: run this code cell under a CUDA-ready environment
%
%Changelog:
%   3/28/2013:  Added processing of History variable if it exists. If
%               multiple files exist for the same data type (i.e. if
%               subject groups were processed in different batches), then
%               the History from the first file will be used.

clc,clear

%add path of GPU function
addpath /nfs/pkg64/gpu/

%user-defined vars---------------------------------------------------------

%general
flagPhase = 1;          %1 = Cue, 2 = probe
flagTTest = 2;          %1 = 1-tailed, 2 = 2-tailed
intRejectLimit = 17;	%maximum number of subjects that can have empty data per electrode
strInputFolder = '/Wavelet results/4-30Hz (full baseline subtraction)/Processed matrices/'; %combined with dirRoot and flagPhase to make the input search path
strOutputFolder = '/Headplots/Permutations/15-30Hz 0-200ms/';    %combined with dirRoot and flagPhase to make save path
boolDoNotOverwrite = true;                                      %useful for multiple-iterations
intRuns = 1;                                                    %number of times to run the same analysis comparisons (repeat mode); note: if boolDoNotOverwrite is not enabled, then this is automatically disabled

%permutation parameters
intCycles = 4000000;  %the number of permutation cycles
dblAlpha = 0.05;    %the alpha value to use in selecting p-values

%frequency settings
arrFrequency = [15 30];

%time scale settings
arrTimes = [0 200];

%comparisons to run
%note: each entry in the structure array is a comparison group with the
%following fieldnames:
%   .name = the name of the comparison, will be used in save filename
%   .pos.dose = the dosage of the first group
%   .pos.cond = the condition of the first group
%   .neg.dose = the dosage of the second group
%   .neg.cond = the condition of the secound group
%*dosages and conditions must match labels used in processed wavelet files
structComparisons(1).name = 'Drug_RvG';
structComparisons(1).pos.dose = 'Drug';
structComparisons(1).pos.cond = 'incong';
structComparisons(1).neg.dose = 'Drug';
structComparisons(1).neg.cond = 'cong';
structComparisons(2).name = 'Plac_RvG';
structComparisons(2).pos.dose = 'Plac';
structComparisons(2).pos.cond = 'incong';
structComparisons(2).neg.dose = 'Plac';
structComparisons(2).neg.cond = 'cong';
structComparisons(3).name = 'Qmoney';
structComparisons(3).pos.dose = 'Drug';
structComparisons(3).pos.cond = 'contrast';
structComparisons(3).neg.dose = 'Plac';
structComparisons(3).neg.cond = 'contrast';

%regions to run
%note: each region has an accompanying switch, turn them on/off to select
%the desired regions to be run
switchRegions = { ...
    'all-Head'       '1';
    };

%subject exclusion
%note: this is used to limit processing to a subset of the subject pool,
%for instance mc-only in the mc-pc dataset
boolUseSubjectExclusion = false;
lstExcludedSubjects = [103,107,110,111,115];

%task-specific vars--------------------------------------------------------

%general
dirRoot = '/nfs/erp-modaf/mc/MC-epoched-Feb11-Glenn/';
strStudy = 'mcpc';                                                      %used in save file name
strMyName = 'pipeHeadplot';                                             %name of this file, used in history
strFilePrefix = 'PRMT-HDPLT';                                           %the prefix added to the beginning of the output filename
pathElectrodeMap = '/nfs/erp-modaf/elec_files/wholehead_elecs.mat';     %defines the regions
intElectrodes = 128;                                                    %max number of electrodes in data

%legacy function call
%note: param1 = row vector of positive term
%      param2 = row vector of negative term
%      param3 = test type
%      param4 = permutation cycles
%      output1 = p-value (range[0,1])
%      output2 = t-value
strFunction = 'permutation_test';   %function name
strTestType = 'paired';             %parameter used to choose permutation method; paired = paired t-test (non-bootstrapping), pooled = mean-difference bootstrapping

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%begin algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Permutation test pipeline: %s\nAuthor: Julian Y. Cheng\n\n',strStudy);

strDataType = 'ERSP';

%parse phase
switch flagPhase
    case 1  %cue
        strPhase = 'Cue';
    case 2  %probe
        strPhase = 'Probe';
    otherwise
        error('incorrect flag set for phase')
end
clear flagPhase

%parse t test tail
switch flagTTest
    case 1  %1-tailed
        boolUseSingleTail = true;
    case 2  %2-tailed
        boolUseSingleTail = false;
    otherwise
        error('incorrect flag set for t test tail')
end
clear flagTTest History

%find processed wavelet data (only load 1st one for parsing parameters, all data loaded after history is saved)
dirInfolder = [dirRoot,strPhase,strInputFolder];
if ~exist(dirInfolder,'dir')
    error('failed to find input folder: %s',dirInfolder)
end
lstInfiles = dir([dirInfolder,strDataType,'*']);
lstInfiles = {lstInfiles.name};
load([dirInfolder,lstInfiles{1}])

%parse frequency and find indexes
arrFrequencyIndexes = find(Data.frequency >= arrFrequency(1),1):find(Data.frequency <= arrFrequency(2),1,'last');

%parse timepoints and find indexes
arrTimesIndexes = find(Data.time >= arrTimes(1),1):find(Data.time <= arrTimes(2),1,'last');

clear Data

%generate output directory (validation after confirm user input)
dirOutput = [dirRoot,strPhase,strOutputFolder];

%construct history
myHistory = fnGetHistory('script',strMyName);

%info messages
fprintf('\n[Parameters]:\n')
fprintf('Input dir: %s\n',dirInfolder)
fprintf('Output dir: %s\n',dirOutput)
fprintf('Phase: %s\n',strPhase)
fprintf('Data type: %s\n',strDataType)
fprintf('Frequency: %i-%iHz\n',arrFrequency(1),arrFrequency(2))
fprintf('Time interval: %i-%ims\n',arrTimes(1),arrTimes(2))
fprintf('Comparisons count: %i\n',length(structComparisons))
fprintf('Regions count: %i\n',sum(str2double(switchRegions(:,2))))
fprintf('Permutation cycles: %i\n',intCycles)
fprintf('Permutation alpha: %.2f\n',dblAlpha)
fprintf('T-test type: %s ',strTestType)
if boolUseSingleTail
    fprintf('1-tailed\n')
else
    fprintf('2-tailed\n')
end
fprintf('Total number of runs: %i\n',intRuns)
fprintf('\n')
if boolUseSubjectExclusion
    fprintf('Warning: subject exclusion is on\n')
end
if ~exist(dirOutput,'dir')
    fprintf('Warning: output folder does not exist and will be automatically created\n')
end
fprintf('\n')

%confirm user-input
strUserInput = input('Continue with the above settings? (Y/N): ','s');
if isempty(strfind(strUserInput,'Y')) && isempty(strfind(strUserInput,'y'))
    error('Operation aborted')
end

%validate output directory
if ~exist(dirOutput,'dir')
    [boolSuccess,strMessage,strMessageID] = mkdir(dirOutput);
    if ~boolSuccess
        error('%s: %s',strMessageID,strMessage)
    end
end

%load in all processed wavelet data MAT files
fprintf('Loading inputs:\n')
for i = 1:length(lstInfiles)
    fprintf('\t%s\n',lstInfiles{i})
    load([dirInfolder,lstInfiles{i}])
    
    %save into data structure
    dataWavelet(i) = Data;
    if exist('History','var')
        dataHistory(i) = History;
    end
    
    clear Data History
end
fprintf('\n')
clear i lstInfiles

%regions definitions-------------------------------------------------------

structRegions(1).name = 'all-Head';
structRegions(1).electrodes = 1:intElectrodes;

%begin processing algorithm------------------------------------------------

fprintf('\n')

%repeat the comparisons structure to run multiple times
for i = 1:intRuns-1
    structComparisons(end+1:end+3) = structComparisons(1:3);
end

%loop through all comparisons
for structComparison = structComparisons
    
    fprintf('Loading data for %s %s Vs %s %s ...', ...
    	structComparison.pos.dose,structComparison.pos.cond,structComparison.neg.dose,structComparison.neg.cond)
    
    %find the data in the processed wavelet data
    idxDosage = find(cellfun(@(x) strcmp(x,structComparison.pos.dose),{dataWavelet.field}));
    idxCondition = find(cellfun(@(x) strcmp(x,structComparison.pos.cond),{dataWavelet.subfield}));
    idxGroup1 = intersect(idxDosage,idxCondition);
    idxDosage = find(cellfun(@(x) strcmp(x,structComparison.neg.dose),{dataWavelet.field}));
    idxCondition = find(cellfun(@(x) strcmp(x,structComparison.neg.cond),{dataWavelet.subfield}));
    idxGroup2 = intersect(idxDosage,idxCondition);
    clearvars idxDosage idxCondition

    %sanity check: make sure subject numbers match
    if ~isequal(dataWavelet(idxGroup1).subject_list,dataWavelet(idxGroup2).subject_list)
        error('subject entries do not match, data is not paired: %s %s Vs %s %s', ...
            structComparison.pos.dose,structComparison.pos.cond,structComparison.neg.dose,structComparison.neg.cond)
    end
    
    %get the data for each group <subject X electrode X frequency X time>
    dataGroup1 = dataWavelet(idxGroup1).matrix;
    dataGroup2 = dataWavelet(idxGroup2).matrix;
    
    %sanity check: make sure both data matrices have the same dimension
    if ~isequal(size(dataGroup1),size(dataGroup2))
        error('data matrices do not match, data is not paired: %s %s Vs %s %s', ...
            structComparison.pos.dose,structComparison.pos.cond,structComparison.neg.dose,structComparison.neg.cond)
    end
    
    %extract data by frequency
    arrFrequencyIndexes = find(dataWavelet(idxGroup1).frequency >= arrFrequency(1),1):find(dataWavelet(idxGroup1).frequency <= arrFrequency(2),1,'last');
    dataGroup1 = dataGroup1(:,:,arrFrequencyIndexes,:);
    arrFrequencyIndexes = find(dataWavelet(idxGroup2).frequency >= arrFrequency(1),1):find(dataWavelet(idxGroup2).frequency <= arrFrequency(2),1,'last');
    dataGroup2 = dataGroup2(:,:,arrFrequencyIndexes,:);
    clear arrFrequencyIndexes

    %extract data by time
    arrTimesIndexes = find(dataWavelet(idxGroup1).time >= arrTimes(1),1):find(dataWavelet(idxGroup1).time <= arrTimes(2),1,'last');
    dataGroup1 = dataGroup1(:,:,:,arrTimesIndexes);
    arrTimesIndexes = find(dataWavelet(idxGroup2).time >= arrTimes(1),1):find(dataWavelet(idxGroup2).time <= arrTimes(2),1,'last');
    dataGroup2 = dataGroup2(:,:,:,arrTimesIndexes);
    clear arrTimesIndexes

    %average down the frequency <subject X electrode X time>
    dataGroup1 = squeeze(nanmean(dataGroup1,3));
    dataGroup2 = squeeze(nanmean(dataGroup2,3));

    %average down the time <subject X electrode>
    dataGroup1 = squeeze(nanmean(dataGroup1,3));
    dataGroup2 = squeeze(nanmean(dataGroup2,3));

    %check if subjection exclusion is enabled
    if boolUseSubjectExclusion
        %find the index of subjects and remove them
        idxDelete = find(ismember(dataWavelet(idxGroup1).subject_list,lstExcludedSubjects));
        dataGroup1(idxDelete,:) = [];
        idxDelete = find(ismember(dataWavelet(idxGroup2).subject_list,lstExcludedSubjects));
        dataGroup2(idxDelete,:) = [];
        
        clearvars idxDelete
    end
    
    fprintf('done\n\n')
    
    %loop through all regions
    for cellRegion = switchRegions' %transpose needed to iterate by rows
        %check if this region should be ran
        if ~str2double(cellRegion{2})
            fprintf('Region:%14s | skipped\n\n',cellRegion{1})
            continue
        end
                
        dataGroup1Subset = dataGroup1;
        dataGroup2Subset = dataGroup2;
        
        %iterate through all frequency and time points
        dataResult.all_t = zeros(1,size(dataGroup1Subset,2));	%used to be called all_t_matrix, in t-values
        dataResult.all_p = zeros(1,size(dataGroup1Subset,2));	%new, in p-values
        dataResult.sig_t = zeros(1,size(dataGroup1Subset,2));	%used to be called sig_matrix, in t-values
        dataResult.sig_p = zeros(1,size(dataGroup1Subset,2));	%new, in p-values
        for i = 1:size(dataGroup1Subset,2)  %electrode
            fprintf('Electrode:%4i | ',i)
            tic
            
            %extract subject data
            arrGroup1 = dataGroup1Subset(:,i);
            arrGroup2 = dataGroup2Subset(:,i);
            
            %find and remove NaN values
            %note: if an electrode is NaN in 1 group, the other group is
            %also excluded because the test requires data to be paired
            lstDeleteIndexes1 = find(isnan(arrGroup1));
            lstDeleteIndexes2 = find(isnan(arrGroup2));
            lstDeleteIndexes = union(lstDeleteIndexes1,lstDeleteIndexes2);
            arrGroup1(lstDeleteIndexes) = [];
            arrGroup2(lstDeleteIndexes) = [];
            clearvars lstDeleteIndexes1 lstDeleteIndexes2
            
            %check if data rejection exceeds threshold
            if (length(lstDeleteIndexes) > intRejectLimit)
                fprintf('Rejected:%3i | EXCLUDED\n',length(lstDeleteIndexes))
                dataResult.all_t(i) = NaN;
                dataResult.all_p(i) = NaN;
                dataResult.sig_t(i) = NaN;
                dataResult.sig_p(i) = NaN;
                [~] = toc;
                continue
            else
                fprintf('Rejected:%3i | Elapsed: ',length(lstDeleteIndexes))
            end
            clearvars lstDeleteIndexes
            
            %call legacy function
            [dblPValue,dblTValue] = eval([strFunction,'(arrGroup1,arrGroup2,strTestType,intCycles)']);
            
            clearvars arrGroup1 arrGroup2
            
            %save the result in all data matrix
            dataResult.all_t(i) = dblTValue;
            dataResult.all_p(i) = dblPValue;
            
            %perform significance comparison
            if boolUseSingleTail
                %assume right tail only
                if (dblPValue > (1-dblAlpha))
                    dataResult.sig_t(i) = dblTValue;
                    dataResult.sig_p(i) = dblPValue;
                end
            else
                if (dblPValue < dblAlpha/2) || (dblPValue > (1-dblAlpha/2))
                    dataResult.sig_t(i) = dblTValue;
                    dataResult.sig_p(i) = dblPValue;
                end
            end
            
            clearvars dblPValue dblTValue
            
            fprintf('%.3f seconds\n',toc)
        end
        
        fprintf('\n')
        
        %create generic variable to store results
        Data = dataResult;
        Data.phase = dataWavelet(idxGroup1).phase;
        Data.type = dataWavelet(idxGroup1).type;
        Data.frequency = dataWavelet(idxGroup1).frequency;
        Data.time = dataWavelet(idxGroup1).time;
        Data.label = structComparison.name;
        Data.region = cellRegion{1};
        clearvars dataResult
        
        %save history information
        if exist('dataHistory','var')
            myHistory.source.Group1_name = dataHistory(idxGroup1).source_name;
            myHistory.source.Group1_id = dataHistory(idxGroup1).source_id;
            myHistory.source.Group2_name = dataHistory(idxGroup2).source_name;
            myHistory.source.Group2_id = dataHistory(idxGroup2).source_id;
        end
        History = myHistory;
        
        %create save filename
        strBuilder = [strFilePrefix,'_'];
        strBuilder = [strBuilder,structComparison.name,'_'];                                    %comparison name
        strBuilder = [strBuilder,num2str(arrFrequency(1)),'-',num2str(arrFrequency(2)),'Hz_'];	%frequency
        strBuilder = [strBuilder,num2str(arrTimes(1)),'-',num2str(arrTimes(2)),'ms'];           %time
        pathOutfile = [dirOutput,strBuilder,'.mat'];
    
        %non-overwrite protection
        intCounter = 2;
        while boolDoNotOverwrite && exist(pathOutfile,'file')
            if isempty(regexp(pathOutfile,'_\d+.mat','match'))
                %add number to end
                pathOutfile = strrep(pathOutfile,'.mat',['_',num2str(intCounter),'.mat']);
                intCounter = intCounter +1;
            else
                %increment the number
                pathOutfile = strrep(pathOutfile,['_',num2str(intCounter-1),'.mat'],['_',num2str(intCounter),'.mat']);
                intCounter = intCounter +1;
            end
        end
        if (intCounter ~= 2)
            fprintf('\nWarning: outfile exists, new file saved as _%i.mat\n',intCounter-1)
        end
        clear intCounter
        
        %save the file
        fprintf('Saving...')
        save(pathOutfile,'Data','History')
        fprintf('done\n\n')
        
        clearvars Data strBuilder dataGroup1Subset dataGroup2Subset pathOutfile
    end
    
    clearvars dataGroup1 dataGroup2
end

clearvars

fprintf('Done\n\n')
beep

%% Single plot (permutation)
%
%This makes plots for all data (all matrices and significant matrices),
%mainly used for debugging. For actual data plotting use multiplot below.

%task-specific vars
pathCleanEEGSet = '/nfs/erp-modaf/mc/ERP/Generic_set_1trial.set';       %needs to have complete chanloc information for all channels
intElectrodes = 128;                                                    %total electrode count, used to construct electrode map of all brain

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%begin algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%request infile
[strFilename,dirFilepath] = uigetfile('','Select permutation MAT file');
if isequal(strFilename,0) || isequal(dirFilepath,   0)
    error('Abort operation')
else
    load([dirFilepath,strFilename])
    
    %data validation
    if ~exist('Data','var') || ~exist('History','var')
        error('Input data is not valid')
    end
end

%get chanlocs from a complete .set file
[~] = evalc('structEEG = pop_loadset(pathCleanEEGSet);');   %silence the function output
structChanlocs = structEEG.chanlocs;

%find electrodes that have been excluded
%note: even though the exclusion list is calculated individually, they
%should be identical
lstElectrodes_master = 1:intElectrodes;
lstIndexes2Delete = find(isnan(Data.all_t));
lstElectrodes_all_t = setdiff(lstElectrodes_master,lstIndexes2Delete);
lstIndexes2Delete = find(isnan(Data.all_p));
lstElectrodes_all_p = setdiff(lstElectrodes_master,lstIndexes2Delete);
lstIndexes2Delete = find(isnan(Data.sig_t));
lstElectrodes_sig_t = setdiff(lstElectrodes_master,lstIndexes2Delete);
lstIndexes2Delete = find(isnan(Data.sig_p));
lstElectrodes_sig_p = setdiff(lstElectrodes_master,lstIndexes2Delete);

%create plot
figure
subplot(2,2,1), colorbar, ...
    h = topoplot(Data.all_t,structChanlocs,'electrodes','on','emarker',{'.','k',5,1},'plotchans',lstElectrodes_all_t);
subplot(2,2,2), colorbar, ...
    h = topoplot(Data.all_p,structChanlocs,'electrodes','on','emarker',{'.','k',5,1},'plotchans',lstElectrodes_all_p);
subplot(2,2,3), colorbar, ...
    h = topoplot(Data.sig_t,structChanlocs,'electrodes','on','emarker',{'.','k',5,1},'plotchans',lstElectrodes_sig_t);
subplot(2,2,4), colorbar, ...
    h = topoplot(Data.sig_p,structChanlocs,'electrodes','on','emarker',{'.','k',5,1},'plotchans',lstElectrodes_sig_p);
    
%% Multi plot (permutation)
%
%This makes all plots in the outfolder under the current directory.
%
%Note: run this code cell after setting the working path to the root of the
%permutation files
%
%Changelog:
%   3/28/2013:  Added processing of History variable if it exists. If
%               multiple files exist for the same data type (i.e. if
%               subject groups were processed in different batches), then
%               the History from the first file will be used.

clear,clc

%user-defined vars
boolPublicationMode = false; %this disables titles and source labeling in the graphs

%task-specific vars
pathCleanEEGSet = '/nfs/erp-modaf/mc/ERP/Generic_set_1trial.set';       %needs to have complete chanloc information for all channels
intElectrodes = 128;                                                    %total electrode count, used to construct electrode map of all brain
strPublicationFolder = 'Publication/';                                  %this is the subfolder that would be created within the normal output directory in publication mode
strFilePrefix = 'PRMT-HDPLT';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%begin algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%find all files to process
lstFiles = dir([strFilePrefix,'*']);
lstFiles = {lstFiles.name};

%check if there are files to process
if isempty(lstFiles)
    error('Cannot find files to process with prefix: %s',strFilePrefix)
end

%create output folders
if boolPublicationMode
    dirOutput_all_t = [strPublicationFolder,'All t figures/'];
    dirOutput_sig_t = [strPublicationFolder,'Significant t figures/'];
    dirOutput_all_p = [strPublicationFolder,'All p figures/'];
    dirOutput_sig_p = [strPublicationFolder,'Significant p figures/'];
    dirOutput_combined = '';    %combined figures for investigational purposes only
    
    if ~exist(dirOutput_all_t,'dir')
        mkdir(dirOutput_all_t);
    end
    if ~exist(dirOutput_sig_t,'dir')
        mkdir(dirOutput_sig_t);
    end
    if ~exist(dirOutput_all_p,'dir')
        mkdir(dirOutput_all_p);
    end
    if ~exist(dirOutput_sig_p,'dir')
        mkdir(dirOutput_sig_p);
    end
else
    dirOutput_all_t = 'All t figures/';
    dirOutput_sig_t = 'Significant t figures/';
    dirOutput_all_p = 'All p figures/';
    dirOutput_sig_p = 'Significant p figures/';
    dirOutput_combined = 'Combined figures/';
    
    if ~exist(dirOutput_all_t,'dir')
        mkdir(dirOutput_all_t);
    end
    if ~exist(dirOutput_sig_t,'dir')
        mkdir(dirOutput_sig_t);
    end
    if ~exist(dirOutput_all_p,'dir')
        mkdir(dirOutput_all_p);
    end
    if ~exist(dirOutput_sig_p,'dir')
        mkdir(dirOutput_sig_p);
    end
    if ~exist(dirOutput_combined,'dir')
        mkdir(dirOutput_combined);
    end
end

%get chanlocs from a complete .set file
[~] = evalc('structEEG = pop_loadset(pathCleanEEGSet);');   %silence the function output
structChanlocs = structEEG.chanlocs;

%process all files
figure
for cellFile = lstFiles
    strFullFilename = cell2mat(cellFile);
    [~,strFilename,~] = fileparts(strFullFilename);
    
    load(strFullFilename)
    
    %data validation
    if ~exist('Data','var') || ~exist('History','var')
        error('Input data not found: %s',strFullFilename)
    end
    
    %note: even though the exclusion list is calculated individually, they
    %should be identical
    lstElectrodes_master = 1:intElectrodes;
    lstIndexes2Delete = find(isnan(Data.all_t));
    lstElectrodes_all_t = setdiff(lstElectrodes_master,lstIndexes2Delete);
    lstIndexes2Delete = find(isnan(Data.all_p));
    lstElectrodes_all_p = setdiff(lstElectrodes_master,lstIndexes2Delete);
    lstIndexes2Delete = find(isnan(Data.sig_t));
    lstElectrodes_sig_t = setdiff(lstElectrodes_master,lstIndexes2Delete);
    lstIndexes2Delete = find(isnan(Data.sig_p));
    lstElectrodes_sig_p = setdiff(lstElectrodes_master,lstIndexes2Delete);
    
    %build source string
    if ismember('source',fieldnames(History))
        strSource1_name = History.source.Group1_name;
        strSource1_id = History.source.Group1_id;
        strSource2_name = History.source.Group1_name;
        strSource2_id = History.source.Group2_id;
        if (strcmp(strSource1_name,strSource2_name) && strcmp(strSource1_id,strSource2_id))
            strSource = ['Source: ',strrep(strSource1_name,'_','\_'),' (',strSource1_id,')'];
        else
            strSource = {['Source: ',strrep(strSource1_name,'_','\_'),' (',strSource1_id,')'],[strrep(strSource2_name,'_','\_'),' (',strSource2_id,')']};
        end
    else
        strSource = ['Source: ',History.vars.dirInfolder];
    end
    
    %Note: not actually plotted in graphs, so commented out to avoid
    %confusion
%     %build title
%     strBuilder = Data.phase;
%     strBuilder = [strBuilder,'-',Data.type];
%     strBuilder = [strBuilder,'-',cell2mat(regexp(strFilename,'\d+-\d+Hz','match'))];
%     strBuilder = [strBuilder,'-',cell2mat(regexp(strFilename,'\d+-\d+ms','match'))];
%     strBuilder = [strBuilder,'-',strrep(Data.label,'_','-')];
    
    %plot individual figures
    
    %all t
    h = topoplot(Data.all_t,structChanlocs,'electrodes','on','emarker',{'.','k',5,1},'plotchans',lstElectrodes_all_t);
    colorbar;
    if ~boolPublicationMode
        lblText = axes('Units','Normal','Position',[0.35 -0.86 .85 .85],'Visible','off');
        set(get(lblText,'Title'),'Visible','on'), title(strSource,'FontSize',6);
    end
    saveas(gcf,[dirOutput_all_t,strFilename,'_2D.tif'],'tiff');
    clf
    
    %sig t
    h = topoplot(Data.sig_t,structChanlocs,'electrodes','on','emarker',{'.','k',5,1},'plotchans',lstElectrodes_sig_t);
    colorbar;
    if ~boolPublicationMode
        lblText = axes('Units','Normal','Position',[0.35 -0.86 .85 .85],'Visible','off');
        set(get(lblText,'Title'),'Visible','on'), title(strSource,'FontSize',6);
    end
    saveas(gcf,[dirOutput_sig_t,strFilename,'_2D.tif'],'tiff');
    clf
    
    %all p
    h = topoplot(Data.all_p,structChanlocs,'electrodes','on','emarker',{'.','k',5,1},'plotchans',lstElectrodes_all_p);
    colorbar;
    if ~boolPublicationMode
        lblText = axes('Units','Normal','Position',[0.35 -0.86 .85 .85],'Visible','off');
        set(get(lblText,'Title'),'Visible','on'), title(strSource,'FontSize',6);
    end
    saveas(gcf,[dirOutput_all_p,strFilename,'_2D.tif'],'tiff');
    clf
    
    %sig p
    h = topoplot(Data.sig_p,structChanlocs,'electrodes','on','emarker',{'.','k',5,1},'plotchans',lstElectrodes_sig_p);
    colorbar;
    if ~boolPublicationMode
        lblText = axes('Units','Normal','Position',[0.35 -0.86 .85 .85],'Visible','off');
        set(get(lblText,'Title'),'Visible','on'), title(strSource,'FontSize',6);
    end
    saveas(gcf,[dirOutput_sig_p,strFilename,'_2D.tif'],'tiff');
    clf

    %combined
    if ~boolPublicationMode
        subplot(2,2,1), colorbar, title('All t'), ...
            h = topoplot(Data.all_t,structChanlocs,'electrodes','on','emarker',{'.','k',5,1},'plotchans',lstElectrodes_all_t);
        subplot(2,2,2), colorbar, title('Sig t'), ...
            h = topoplot(Data.sig_t,structChanlocs,'electrodes','on','emarker',{'.','k',5,1},'plotchans',lstElectrodes_sig_t);
        subplot(2,2,3), colorbar, title('All p'), ...
            h = topoplot(Data.all_p,structChanlocs,'electrodes','on','emarker',{'.','k',5,1},'plotchans',lstElectrodes_all_p);
        subplot(2,2,4), colorbar, title('Sig p'), ...
            h = topoplot(Data.sig_p,structChanlocs,'electrodes','on','emarker',{'.','k',5,1},'plotchans',lstElectrodes_sig_p);
            lblText = axes('Units','Normal','Position',[0.35 -0.86 .85 .85],'Visible','off');
            set(get(lblText,'Title'),'Visible','on'), title(strSource,'FontSize',6);
        saveas(gcf,['Combined figures/',strFilename,'_2D.tif'],'tiff');
        clf
    end
    
    clear Data History
end

close,beep