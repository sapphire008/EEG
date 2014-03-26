%Spectrogram plot generation sub-routine
%Author: Julian Y. Cheng
%2/19/2013
%
%This script is an adaption of the Create_Spectrograms_modaf.m in order to
%make use of the new condition MAT files. Not all functionalities are
%implemented; namely, the electrode-by-electrode graphing is not
%implemented because it was for investigational purposes only.

%% Single plot
%
%This code cell processes a single input file per user selection. Used for
%debugging purposes. First attempt to implement UI-based switch selection
%
%Changelog:
%   3/21/2013:  Added processing of early VS late analysis. Note that for
%               new code the syntax has been changed to be more generic
%               (ex. Dosage -> field); this has been done after realizing
%               that other datasets may not have drug/condition combination
%   3/25/2013:  Added processing of block-wise analysis.

clc,clear

%user-defined vars

%color scale
boolUseCustomScaling = true;                        %for automatic-scaling, set to false
structColorScale(1).type = 'ERSP';                  %if type is changed, code modification is necessary
structColorScale(1).limits = [-0.33 0.33];
structColorScale(2).type = 'ITC single condition';  %if type is changed, code modification is necessary
structColorScale(2).limits = [0 0.3];
structColorScale(3).type = 'ITC contrasts';         %if type is changed, code modification is necessary
structColorScale(3).limits = [-0.04 0.04];

%subject exclusion
boolUseSubjectExclusion = false;
lstExcludedSubjects = [103,108,110,111,115];   %will not process these subjects, use to limit the MC/PC dataset to MC or PC only

%task-specific vars

%general
dirRoot = '/nfs/erp-modaf/mc/MC-epoched-Feb11-Glenn/';
strInputFolder = 'Wavelet results/';                                    %combined with dirRoot and phase string to produce path to input folders
strSubfolder = 'Processed matrices/';                                   %the folder to search for the input MAT files
pathElectrodeMap = '/nfs/erp-modaf/elec_files/wholehead_elecs.mat';     %defines the regions
lstPhase = {'Cue','Probe'};


%comparisons
%note: this specifies which comparisons are graphable, by defining the
%required files here. Note that for Qmoney conditions, the order of the
%files define the order of subtraction; i.e. file1 - file2. If 4 files are
%defined, then the algorithm becomes: (file1 - file2) - (file3 - file4).
%note: the "if true...end" statements are used for code folding.
%[all blocks]
if true
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
end
%[early blocks]
if true
    structComparisons(10).name = 'Drug green (early)';
    structComparisons(10).files = {'Drug-early_cong'}; 
    structComparisons(11).name = 'Drug red (early)';
    structComparisons(11).files = {'Drug-early_incong'}; 
    structComparisons(12).name = 'Drug contrast (early)';
    structComparisons(12).files = {'Drug-early_contrast'}; 
    structComparisons(13).name = 'Plac green (early)';
    structComparisons(13).files = {'Plac-early_cong'}; 
    structComparisons(14).name = 'Plac red (early)';
    structComparisons(14).files = {'Plac-early_incong'}; 
    structComparisons(15).name = 'Plac contrast (early)';
    structComparisons(15).files = {'Plac-early_contrast'}; 
    structComparisons(16).name = 'Drug green - Plac green (early)';
    structComparisons(16).files = {'Drug-early_cong','Plac-early_cong'}; 
    structComparisons(17).name = 'Drug red - Plac red (early)';
    structComparisons(17).files = {'Drug-early_incong','Plac-early_incong'}; 
    structComparisons(18).name = 'Qmoney (early)';
    structComparisons(18).files = {'Drug-early_contrast','Plac-early_contrast'}; 
end
%[late blocks]
if true
    structComparisons(19).name = 'Drug green (late)';
    structComparisons(19).files = {'Drug-late_cong'}; 
    structComparisons(20).name = 'Drug red (late)';
    structComparisons(20).files = {'Drug-late_incong'}; 
    structComparisons(21).name = 'Drug contrast (late)';
    structComparisons(21).files = {'Drug-late_contrast'}; 
    structComparisons(22).name = 'Plac green (late)';
    structComparisons(22).files = {'Plac-late_cong'}; 
    structComparisons(23).name = 'Plac red (late)';
    structComparisons(23).files = {'Plac-late_incong'}; 
    structComparisons(24).name = 'Plac contrast (late)';
    structComparisons(24).files = {'Plac-late_contrast'}; 
    structComparisons(25).name = 'Drug green - Plac green (late)';
    structComparisons(25).files = {'Drug-late_cong','Plac-late_cong'}; 
    structComparisons(26).name = 'Drug red - Plac red (late)';
    structComparisons(26).files = {'Drug-late_incong','Plac-late_incong'}; 
    structComparisons(27).name = 'Qmoney (late)';
    structComparisons(27).files = {'Drug-late_contrast','Plac-late_contrast'}; 
end
%[late VS early]
if true
    structComparisons(28).name = 'Drug green (late-VS-early)';
    structComparisons(28).files = {'Drug-late_cong','Drug-early_cong'}; 
    structComparisons(29).name = 'Drug red (late-VS-early)';
    structComparisons(29).files = {'Drug-late_incong','Drug-early_incong'}; 
    structComparisons(30).name = 'Drug contrast (late-VS-early)';
    structComparisons(30).files = {'Drug-late_contrast','Drug-early_contrast'}; 
    structComparisons(31).name = 'Plac green (late-VS-early)';
    structComparisons(31).files = {'Plac-late_cong','Plac-early_cong'}; 
    structComparisons(32).name = 'Plac red (late-VS-early)';
    structComparisons(32).files = {'Plac-late_incong','Plac-early_incong'}; 
    structComparisons(33).name = 'Plac contrast (late-VS-early)';
    structComparisons(33).files = {'Plac-late_contrast','Plac-early_contrast'}; 
    structComparisons(34).name = 'Drug green - Plac green (late-VS-early)';
    structComparisons(34).files = {'Drug-late_cong','Plac-late_cong','Drug-early_cong','Plac-early_cong'}; 
    structComparisons(35).name = 'Drug red - Plac red (late-VS-early)';
    structComparisons(35).files = {'Drug-late_incong','Plac-late_incong','Drug-early_incong','Plac-early_incong'}; 
    structComparisons(36).name = 'Qmoney (late-VS-early)';
    structComparisons(36).files = {'Drug-late_contrast','Plac-late_contrast','Drug-early_contrast','Plac-early_contrast'};
end
%[blocks 1,2]
if true
    structComparisons(37).name = 'Drug green (blocks 12)';
    structComparisons(37).files = {'Drug-blk12_cong'}; 
    structComparisons(38).name = 'Drug red (blocks 12)';
    structComparisons(38).files = {'Drug-blk12_incong'}; 
    structComparisons(39).name = 'Drug contrast (blocks 12)';
    structComparisons(39).files = {'Drug-blk12_contrast'}; 
    structComparisons(40).name = 'Plac green (blocks 12)';
    structComparisons(40).files = {'Plac-blk12_cong'}; 
    structComparisons(41).name = 'Plac red (blocks 12)';
    structComparisons(41).files = {'Plac-blk12_incong'}; 
    structComparisons(42).name = 'Plac contrast (blocks 12)';
    structComparisons(42).files = {'Plac-blk12_contrast'}; 
    structComparisons(43).name = 'Drug green - Plac green (blocks 12)';
    structComparisons(43).files = {'Drug-blk12_cong','Plac-blk12_cong'}; 
    structComparisons(44).name = 'Drug red - Plac red (blocks 12)';
    structComparisons(44).files = {'Drug-blk12_incong','Plac-blk12_incong'}; 
    structComparisons(45).name = 'Qmoney (blocks 12)';
    structComparisons(45).files = {'Drug-blk12_contrast','Plac-blk12_contrast'}; 
end
%[blocks 3,4]
if true
    structComparisons(46).name = 'Drug green (blocks 34)';
    structComparisons(46).files = {'Drug-blk34_cong'}; 
    structComparisons(47).name = 'Drug red (blocks 34)';
    structComparisons(47).files = {'Drug-blk34_incong'}; 
    structComparisons(48).name = 'Drug contrast (blocks 34)';
    structComparisons(48).files = {'Drug-blk34_contrast'}; 
    structComparisons(49).name = 'Plac green (blocks 34)';
    structComparisons(49).files = {'Plac-blk34_cong'}; 
    structComparisons(50).name = 'Plac red (blocks 34)';
    structComparisons(50).files = {'Plac-blk34_incong'}; 
    structComparisons(51).name = 'Plac contrast (blocks 34)';
    structComparisons(51).files = {'Plac-blk34_contrast'}; 
    structComparisons(52).name = 'Drug green - Plac green (blocks 34)';
    structComparisons(52).files = {'Drug-blk34_cong','Plac-blk34_cong'}; 
    structComparisons(53).name = 'Drug red - Plac red (blocks 34)';
    structComparisons(53).files = {'Drug-blk34_incong','Plac-blk34_incong'}; 
    structComparisons(54).name = 'Qmoney (blocks 34)';
    structComparisons(54).files = {'Drug-blk34_contrast','Plac-blk34_contrast'}; 
end
%[blocks 5,6]
if true
    structComparisons(55).name = 'Drug green (blocks 56)';
    structComparisons(55).files = {'Drug-blk56_cong'}; 
    structComparisons(56).name = 'Drug red (blocks 56)';
    structComparisons(56).files = {'Drug-blk56_incong'}; 
    structComparisons(57).name = 'Drug contrast (blocks 56)';
    structComparisons(57).files = {'Drug-blk56_contrast'}; 
    structComparisons(58).name = 'Plac green (blocks 56)';
    structComparisons(58).files = {'Plac-blk56_cong'}; 
    structComparisons(59).name = 'Plac red (blocks 56)';
    structComparisons(59).files = {'Plac-blk56_incong'}; 
    structComparisons(60).name = 'Plac contrast (blocks 56)';
    structComparisons(60).files = {'Plac-blk56_contrast'}; 
    structComparisons(61).name = 'Drug green - Plac green (blocks 56)';
    structComparisons(61).files = {'Drug-blk56_cong','Plac-blk56_cong'}; 
    structComparisons(62).name = 'Drug red - Plac red (blocks 56)';
    structComparisons(62).files = {'Drug-blk56_incong','Plac-blk56_incong'}; 
    structComparisons(63).name = 'Qmoney (blocks 56)';
    structComparisons(63).files = {'Drug-blk56_contrast','Plac-blk56_contrast'}; 
end
%[blocks 7,8]
if true
    structComparisons(64).name = 'Drug green (blocks 78)';
    structComparisons(64).files = {'Drug-blk78_cong'}; 
    structComparisons(65).name = 'Drug red (blocks 78)';
    structComparisons(65).files = {'Drug-blk78_incong'}; 
    structComparisons(66).name = 'Drug contrast (blocks 78)';
    structComparisons(66).files = {'Drug-blk78_contrast'}; 
    structComparisons(67).name = 'Plac green (blocks 78)';
    structComparisons(67).files = {'Plac-blk78_cong'}; 
    structComparisons(68).name = 'Plac red (blocks 78)';
    structComparisons(68).files = {'Plac-blk78_incong'}; 
    structComparisons(69).name = 'Plac contrast (blocks 78)';
    structComparisons(69).files = {'Plac-blk78_contrast'}; 
    structComparisons(70).name = 'Drug green - Plac green (blocks 78)';
    structComparisons(70).files = {'Drug-blk78_cong','Plac-blk78_cong'}; 
    structComparisons(71).name = 'Drug red - Plac red (blocks 78)';
    structComparisons(71).files = {'Drug-blk78_incong','Plac-blk78_incong'}; 
    structComparisons(72).name = 'Qmoney (blocks 78)';
    structComparisons(72).files = {'Drug-blk78_contrast','Plac-blk78_contrast'}; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%begin algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%regions definitions-------------------------------------------------------

load(pathElectrodeMap)

structRegions(1).name = 'all-PFC';
structRegions(1).electrodes = cat(2, wholehead_elecs{1,2},wholehead_elecs{2,2}, wholehead_elecs{3,2});
structRegions(2).name = 'left-PFC';
structRegions(2).electrodes = wholehead_elecs{2,2};
structRegions(3).name = 'mid-PFC';
structRegions(3).electrodes = wholehead_elecs{1,2};
structRegions(4).name = 'right-PFC';
structRegions(4).electrodes = wholehead_elecs{3,2};
structRegions(5).name = 'all-Central';
structRegions(5).electrodes = cat(2, wholehead_elecs{8,2},wholehead_elecs{4,2}, wholehead_elecs{9,2});
structRegions(6).name = 'left-Central';
structRegions(6).electrodes = wholehead_elecs{8,2};
structRegions(7).name = 'mid-Central';
structRegions(7).electrodes = wholehead_elecs{4,2};
structRegions(8).name = 'right-Central';
structRegions(8).electrodes = wholehead_elecs{9,2};
structRegions(9).name = 'all-Post';
structRegions(9).electrodes = cat(2, wholehead_elecs{6,2},wholehead_elecs{5,2}, wholehead_elecs{7,2});
structRegions(10).name = 'left-Post';
structRegions(10).electrodes = wholehead_elecs{6,2};
structRegions(11).name = 'mid-Post';
structRegions(11).electrodes = wholehead_elecs{5,2};
structRegions(12).name = 'right-Post';
structRegions(12).electrodes = wholehead_elecs{7,2};

clear wholehead_elecs

%begin ui algorithm--------------------------------------------------------

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

%find ERSP/ITC datatypes
boolFoundERSP = any(cell2mat(strfind(lstInfiles,'ERSP')));
boolFoundITC = any(cell2mat(strfind(lstInfiles,'ITC')));

%create selection list
lstDatatype = {};
if boolFoundERSP
    lstDatatype{end+1} = 'ERSP';
end
if boolFoundITC
    lstDatatype{end+1} = 'ITC';
end
if boolFoundERSP && boolFoundITC
    lstDatatype{end+1} = 'ERSP & ITC';
end
clear boolFoundERSP boolFoundITC

%show dialog
[idxSelection boolOK] = listdlg('Name','Data type','PromptString','Select data type(s):', ...
                                'ListString',lstDatatype,'ListSize',[200 300],'SelectionMode','single');
if boolOK
    switch lstDatatype{idxSelection}
        case 'ERSP'
            idxDelete = find(~cellfun(@isempty,strfind(lstInfiles,'ITC')));
            lstInfiles(idxDelete) = [];
        case 'ITC'
            idxDelete = find(~cellfun(@isempty,strfind(lstInfiles,'ERSP')));
            lstInfiles(idxDelete) = [];
        case 'ERSP & ITC'
            %do nothing
        otherwise
            error('Failed to parse selection for datatype; undefined selection: %s',lstDatatype{idxSelection})
    end
else
    error('Operation aborted')
end
clear lstDatatype idxSelection boolOK

%prompt for comparisons to graph

%cleanup comparisons to only the ones that have required files found
lstDeleteIndexes = [];
for i = 1:length(structComparisons)
    for j = 1:length(structComparisons(i).files)
        if all(cellfun(@isempty,regexp(lstInfiles,structComparisons(i).files(j))))
            %exclude this comparison if any required file is not found
            lstDeleteIndexes(end+1) = i;
        end
    end
end
structComparisons(lstDeleteIndexes) = [];
clear lstDeleteIndexes i j

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
    %Update 3/21/2013: the numbers will double further if performing one of
    %the 4 file subtractions for late VS early analysis
    lstFileIndexes = cell2mat(lstFileIndexes(:)');  %discard sorting by search string, instead just load all files like in pipePermutationTest
    for i = 1:length(lstFileIndexes)
        load([dirAnalysis,lstInfiles{lstFileIndexes(i)}])
        
        %save into data structure
        dataWavelet(i) = Data;
        
        clear Data
    end
    clear lstFileIndexes
    
    %loop through all regions
    for myRegion = structRegions
        %get the electrode list
        lstElectrodes = myRegion.electrodes;
        
        %process both ERSP and ITC
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
            if (length(idxDatatype) == 4)   %(A-B)-(C-D)
                %parse dosage and condition information from myComparison.files
                strSplit = regexp(myComparison.files{1},'_','split');
                strField1 = strSplit{1};
                strSubfield1 = strSplit{2};
                strSplit = regexp(myComparison.files{2},'_','split');
                strField2 = strSplit{1};
                strSubfield2 = strSplit{2};
                strSplit = regexp(myComparison.files{3},'_','split');
                strField3 = strSplit{1};
                strSubfield3 = strSplit{2};
                strSplit = regexp(myComparison.files{4},'_','split');
                strField4 = strSplit{1};
                strSubfield4 = strSplit{2};
                
                %find where the conditions are stored in the dataWavelet structure
                idxField = find(cellfun(@(x) strcmp(x,strField1),{dataWavelet.field}));
                idxSubfield = find(cellfun(@(x) strcmp(x,strSubfield1),{dataWavelet.subfield}));
                idxGroup1 = intersect(idxField,idxSubfield);
                idxGroup1 = intersect(idxDatatype,idxGroup1);
                idxField = find(cellfun(@(x) strcmp(x,strField2),{dataWavelet.field}));
                idxSubfield = find(cellfun(@(x) strcmp(x,strSubfield2),{dataWavelet.subfield}));
                idxGroup2 = intersect(idxField,idxSubfield);
                idxGroup2 = intersect(idxDatatype,idxGroup2);
                idxField = find(cellfun(@(x) strcmp(x,strField3),{dataWavelet.field}));
                idxSubfield = find(cellfun(@(x) strcmp(x,strSubfield3),{dataWavelet.subfield}));
                idxGroup3 = intersect(idxField,idxSubfield);
                idxGroup3 = intersect(idxDatatype,idxGroup3);
                idxField = find(cellfun(@(x) strcmp(x,strField4),{dataWavelet.field}));
                idxSubfield = find(cellfun(@(x) strcmp(x,strSubfield4),{dataWavelet.subfield}));
                idxGroup4 = intersect(idxField,idxSubfield);
                idxGroup4 = intersect(idxDatatype,idxGroup4);
                
                %perform subtraction, and modify structure of Data
                myData = dataWavelet(idxGroup1);
                myData.field = myComparison.name;
                myData.subfield = '';   %this is ignored
                myData.matrix = (dataWavelet(idxGroup1).matrix - dataWavelet(idxGroup2).matrix) - ...
                                (dataWavelet(idxGroup3).matrix - dataWavelet(idxGroup4).matrix);
                
                clear strField1 strField2 strField3 strField4 strSubfield1 strSubfield2 strSubfield3 strSubfield4
                clear idxField idxSubfield idxGroup1 idxGroup2 idxGroup3 idxGroup4
            elseif (length(idxDatatype) == 2)   %(A-B)
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
                
                clear strDosage1 strDosage2 strCondition1 strCondition2
                clear idxDosage idxCondition idxGroup1 idxGroup2
            elseif (length(idxDatatype) == 1)   %(A-Null)
                myData = dataWavelet(idxDatatype);
            else
                error('Invalid number of datasets found (should be 1,2,4): %i',length(idxDatatype))
            end
                
            %check if subject exclusion is enabled
            if boolUseSubjectExclusion
                idxDelete = find(ismember(myData.subject_list,lstExcludedSubjects));
                myData.matrix(idxDelete,:,:,:) = [];
                
                clear idxDelete
            end
            
            %extract electrode data
            myData.matrix = myData.matrix(:,lstElectrodes,:,:);
            
            %average down the region to 1 giant electrode <subject X frequency X time>
            myData.matrix = squeeze(nanmean(myData.matrix,2));
            
            %average down all subjects to get group level results <frequency X time>
            myData.matrix = squeeze(nanmean(myData.matrix,1));
            
            %set color scale
            if strcmp(strDatatype,'ERSP')
                idxColorScale = find(cellfun(@(x) strcmp(x,'ERSP'),{structColorScale.type}));
            elseif (min(myData.matrix(:)) >= 0)
                idxColorScale = find(cellfun(@(x) strcmp(x,'ITC single condition'),{structColorScale.type}));
            else
                idxColorScale = find(cellfun(@(x) strcmp(x,'ITC contrasts'),{structColorScale.type}));
            end
            
            %graph the result
            figure
            h = tftopo(myData.matrix,myData.time,myData.frequency,'title',['{\bf{\color{blue}',strDatatype,'}} ',myComparison.name,' {\color{red}',myRegion.name,'}'],'verbose','off');
            if boolUseCustomScaling
                caxis(structColorScale(idxColorScale).limits), ...
                colorbar;
            else
                colorbar;
            end
            lblText = axes('Units','Normal','Position',[0.3 -0.86 .85 .85],'Visible','off');
            set(get(lblText,'Title'),'Visible','on'), title(['Source: ',strrep(myData.source,'_','\_')],'FontSize',8);
            
            clear myData idxDatatype strDatatype lblText idxColorScale 
        end
        clear lstElectrodes
    end
end

%% Multi plot
%
%This code cell is used to graph all plots from a single input folder. Uses
%UI based switching for phase, analysis folder selection, regions, and 
%output folder specification. First implementation of textboxes.
%
%Changelog:
%   3/21/2013:  Added processing of early VS late analysis. Note that for
%               new code the syntax has been changed to be more generic
%               (ex. Dosage -> field); this has been done after realizing
%               that other datasets may not have drug/condition combination
%   3/25/2013:  Added processing of block-wise analysis.
%   3/28/2013:  Added processing of History variable if it exists. If
%               multiple files exist for the same data type (i.e. if
%               subject groups were processed in different batches), then
%               the History from the first file will be used.

clc,close,clear

%user-defined vars

boolPublicationMode = false; %this disables titles and source labeling in the graphs

%color scale
boolUseCustomScaling = true;                        %for automatic-scaling, set to false
structColorScale(1).type = 'ERSP';                  %if type is changed, code modification is necessary
structColorScale(1).limits = [-0.33 0.33];
structColorScale(2).type = 'ITC single condition';  %if type is changed, code modification is necessary
structColorScale(2).limits = [0 0.3];
structColorScale(3).type = 'ITC contrasts';         %if type is changed, code modification is necessary
structColorScale(3).limits = [-0.04 0.04];

%subject exclusion
boolUseSubjectExclusion = false;
lstExcludedSubjects = [103,108,110,111,115];   %will not process these subjects, use to limit the MC/PC dataset to MC or PC only

%task-specific vars

%general
dirRoot = '/nfs/erp-modaf/ms/MS-epoched-Feb11-Glenn/';
strInputFolder = 'Wavelet results/';                                    %combined with dirRoot and phase string to produce path to input folders
strSubfolder = 'Processed matrices/';                                   %the folder to search for the input MAT files
strOutputFolder = 'Spectrograms/';
strPublicationFolder = 'Publication/';                                  %this is the subfolder that would be created within the normal output directory in publication mode
strFilePrefix = 'SPCT';
pathElectrodeMap = '/nfs/erp-modaf/elec_files/wholehead_elecs.mat';     %defines the regions
lstPhase = {'Cue','Probe'};

%comparisons
%note: this specifies which comparisons are graphable, by defining the
%required files here. Note that for Qmoney conditions, the order of the
%files define the order of subtraction; i.e. file1 - file2. If 4 files are
%defined, then the algorithm becomes: (file1 - file2) - (file3 - file4).
%note: the "if true...end" statements are used for code folding
%[all blocks]
if true
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
end
%[early blocks]
if true
    structComparisons(10).name = 'Drug green (early)';
    structComparisons(10).files = {'Drug-early_cong'}; 
    structComparisons(11).name = 'Drug red (early)';
    structComparisons(11).files = {'Drug-early_incong'}; 
    structComparisons(12).name = 'Drug contrast (early)';
    structComparisons(12).files = {'Drug-early_contrast'}; 
    structComparisons(13).name = 'Plac green (early)';
    structComparisons(13).files = {'Plac-early_cong'}; 
    structComparisons(14).name = 'Plac red (early)';
    structComparisons(14).files = {'Plac-early_incong'}; 
    structComparisons(15).name = 'Plac contrast (early)';
    structComparisons(15).files = {'Plac-early_contrast'}; 
    structComparisons(16).name = 'Drug green - Plac green (early)';
    structComparisons(16).files = {'Drug-early_cong','Plac-early_cong'}; 
    structComparisons(17).name = 'Drug red - Plac red (early)';
    structComparisons(17).files = {'Drug-early_incong','Plac-early_incong'}; 
    structComparisons(18).name = 'Qmoney (early)';
    structComparisons(18).files = {'Drug-early_contrast','Plac-early_contrast'}; 
end
%[late blocks]
if true
    structComparisons(19).name = 'Drug green (late)';
    structComparisons(19).files = {'Drug-late_cong'}; 
    structComparisons(20).name = 'Drug red (late)';
    structComparisons(20).files = {'Drug-late_incong'}; 
    structComparisons(21).name = 'Drug contrast (late)';
    structComparisons(21).files = {'Drug-late_contrast'}; 
    structComparisons(22).name = 'Plac green (late)';
    structComparisons(22).files = {'Plac-late_cong'}; 
    structComparisons(23).name = 'Plac red (late)';
    structComparisons(23).files = {'Plac-late_incong'}; 
    structComparisons(24).name = 'Plac contrast (late)';
    structComparisons(24).files = {'Plac-late_contrast'}; 
    structComparisons(25).name = 'Drug green - Plac green (late)';
    structComparisons(25).files = {'Drug-late_cong','Plac-late_cong'}; 
    structComparisons(26).name = 'Drug red - Plac red (late)';
    structComparisons(26).files = {'Drug-late_incong','Plac-late_incong'}; 
    structComparisons(27).name = 'Qmoney (late)';
    structComparisons(27).files = {'Drug-late_contrast','Plac-late_contrast'}; 
end
%[late VS early]
if true
    structComparisons(28).name = 'Drug green (late-VS-early)';
    structComparisons(28).files = {'Drug-late_cong','Drug-early_cong'}; 
    structComparisons(29).name = 'Drug red (late-VS-early)';
    structComparisons(29).files = {'Drug-late_incong','Drug-early_incong'}; 
    structComparisons(30).name = 'Drug contrast (late-VS-early)';
    structComparisons(30).files = {'Drug-late_contrast','Drug-early_contrast'}; 
    structComparisons(31).name = 'Plac green (late-VS-early)';
    structComparisons(31).files = {'Plac-late_cong','Plac-early_cong'}; 
    structComparisons(32).name = 'Plac red (late-VS-early)';
    structComparisons(32).files = {'Plac-late_incong','Plac-early_incong'}; 
    structComparisons(33).name = 'Plac contrast (late-VS-early)';
    structComparisons(33).files = {'Plac-late_contrast','Plac-early_contrast'}; 
    structComparisons(34).name = 'Drug green - Plac green (late-VS-early)';
    structComparisons(34).files = {'Drug-late_cong','Plac-late_cong','Drug-early_cong','Plac-early_cong'}; 
    structComparisons(35).name = 'Drug red - Plac red (late-VS-early)';
    structComparisons(35).files = {'Drug-late_incong','Plac-late_incong','Drug-early_incong','Plac-early_incong'}; 
    structComparisons(36).name = 'Qmoney (late-VS-early)';
    structComparisons(36).files = {'Drug-late_contrast','Plac-late_contrast','Drug-early_contrast','Plac-early_contrast'};
end
%[blocks 1,2]
if true
    structComparisons(37).name = 'Drug green (blocks 12)';
    structComparisons(37).files = {'Drug-blk12_cong'}; 
    structComparisons(38).name = 'Drug red (blocks 12)';
    structComparisons(38).files = {'Drug-blk12_incong'}; 
    structComparisons(39).name = 'Drug contrast (blocks 12)';
    structComparisons(39).files = {'Drug-blk12_contrast'}; 
    structComparisons(40).name = 'Plac green (blocks 12)';
    structComparisons(40).files = {'Plac-blk12_cong'}; 
    structComparisons(41).name = 'Plac red (blocks 12)';
    structComparisons(41).files = {'Plac-blk12_incong'}; 
    structComparisons(42).name = 'Plac contrast (blocks 12)';
    structComparisons(42).files = {'Plac-blk12_contrast'}; 
    structComparisons(43).name = 'Drug green - Plac green (blocks 12)';
    structComparisons(43).files = {'Drug-blk12_cong','Plac-blk12_cong'}; 
    structComparisons(44).name = 'Drug red - Plac red (blocks 12)';
    structComparisons(44).files = {'Drug-blk12_incong','Plac-blk12_incong'}; 
    structComparisons(45).name = 'Qmoney (blocks 12)';
    structComparisons(45).files = {'Drug-blk12_contrast','Plac-blk12_contrast'}; 
end
%[blocks 3,4]
if true
    structComparisons(46).name = 'Drug green (blocks 34)';
    structComparisons(46).files = {'Drug-blk34_cong'}; 
    structComparisons(47).name = 'Drug red (blocks 34)';
    structComparisons(47).files = {'Drug-blk34_incong'}; 
    structComparisons(48).name = 'Drug contrast (blocks 34)';
    structComparisons(48).files = {'Drug-blk34_contrast'}; 
    structComparisons(49).name = 'Plac green (blocks 34)';
    structComparisons(49).files = {'Plac-blk34_cong'}; 
    structComparisons(50).name = 'Plac red (blocks 34)';
    structComparisons(50).files = {'Plac-blk34_incong'}; 
    structComparisons(51).name = 'Plac contrast (blocks 34)';
    structComparisons(51).files = {'Plac-blk34_contrast'}; 
    structComparisons(52).name = 'Drug green - Plac green (blocks 34)';
    structComparisons(52).files = {'Drug-blk34_cong','Plac-blk34_cong'}; 
    structComparisons(53).name = 'Drug red - Plac red (blocks 34)';
    structComparisons(53).files = {'Drug-blk34_incong','Plac-blk34_incong'}; 
    structComparisons(54).name = 'Qmoney (blocks 34)';
    structComparisons(54).files = {'Drug-blk34_contrast','Plac-blk34_contrast'}; 
end
%[blocks 5,6]
if true
    structComparisons(55).name = 'Drug green (blocks 56)';
    structComparisons(55).files = {'Drug-blk56_cong'}; 
    structComparisons(56).name = 'Drug red (blocks 56)';
    structComparisons(56).files = {'Drug-blk56_incong'}; 
    structComparisons(57).name = 'Drug contrast (blocks 56)';
    structComparisons(57).files = {'Drug-blk56_contrast'}; 
    structComparisons(58).name = 'Plac green (blocks 56)';
    structComparisons(58).files = {'Plac-blk56_cong'}; 
    structComparisons(59).name = 'Plac red (blocks 56)';
    structComparisons(59).files = {'Plac-blk56_incong'}; 
    structComparisons(60).name = 'Plac contrast (blocks 56)';
    structComparisons(60).files = {'Plac-blk56_contrast'}; 
    structComparisons(61).name = 'Drug green - Plac green (blocks 56)';
    structComparisons(61).files = {'Drug-blk56_cong','Plac-blk56_cong'}; 
    structComparisons(62).name = 'Drug red - Plac red (blocks 56)';
    structComparisons(62).files = {'Drug-blk56_incong','Plac-blk56_incong'}; 
    structComparisons(63).name = 'Qmoney (blocks 56)';
    structComparisons(63).files = {'Drug-blk56_contrast','Plac-blk56_contrast'}; 
end
%[blocks 7,8]
if true
    structComparisons(64).name = 'Drug green (blocks 78)';
    structComparisons(64).files = {'Drug-blk78_cong'}; 
    structComparisons(65).name = 'Drug red (blocks 78)';
    structComparisons(65).files = {'Drug-blk78_incong'}; 
    structComparisons(66).name = 'Drug contrast (blocks 78)';
    structComparisons(66).files = {'Drug-blk78_contrast'}; 
    structComparisons(67).name = 'Plac green (blocks 78)';
    structComparisons(67).files = {'Plac-blk78_cong'}; 
    structComparisons(68).name = 'Plac red (blocks 78)';
    structComparisons(68).files = {'Plac-blk78_incong'}; 
    structComparisons(69).name = 'Plac contrast (blocks 78)';
    structComparisons(69).files = {'Plac-blk78_contrast'}; 
    structComparisons(70).name = 'Drug green - Plac green (blocks 78)';
    structComparisons(70).files = {'Drug-blk78_cong','Plac-blk78_cong'}; 
    structComparisons(71).name = 'Drug red - Plac red (blocks 78)';
    structComparisons(71).files = {'Drug-blk78_incong','Plac-blk78_incong'}; 
    structComparisons(72).name = 'Qmoney (blocks 78)';
    structComparisons(72).files = {'Drug-blk78_contrast','Plac-blk78_contrast'}; 
end
%[MS dose 3]
%note: causes issues because not 1-1, need to rethink this
if false
    structComparisons(73).name = 'Drug green (dose 3)';
    structComparisons(73).files = {'Drug-dose3_cong'}; 
    structComparisons(74).name = 'Drug red (dose 3)';
    structComparisons(74).files = {'Drug-dose3_incong'}; 
    structComparisons(75).name = 'Drug contrast (dose 3)';
    structComparisons(75).files = {'Drug-dose3_contrast'}; 
    structComparisons(76).name = 'Plac green (dose 3)';
    structComparisons(76).files = {'Plac-dose3_cong'}; 
    structComparisons(77).name = 'Plac red (dose 3)';
    structComparisons(77).files = {'Plac-dose3_incong'}; 
    structComparisons(78).name = 'Plac contrast (dose 3)';
    structComparisons(78).files = {'Plac-dose3_contrast'}; 
    structComparisons(79).name = 'Drug green - Plac green (dose 3)';
    structComparisons(79).files = {'Drug-dose3_cong','Plac-dose3_cong'}; 
    structComparisons(80).name = 'Drug red - Plac red (dose 3)';
    structComparisons(80).files = {'Drug-dose3_incong','Plac-dose3_incong'}; 
    structComparisons(81).name = 'Qmoney (dose 3)';
    structComparisons(81).files = {'Drug-dose3_contrast','Plac-dose3_contrast'}; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%begin algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%regions definitions-------------------------------------------------------

load(pathElectrodeMap)

structRegions(1).name = 'all-PFC';
structRegions(1).electrodes = cat(2, wholehead_elecs{1,2},wholehead_elecs{2,2}, wholehead_elecs{3,2});
structRegions(2).name = 'left-PFC';
structRegions(2).electrodes = wholehead_elecs{2,2};
structRegions(3).name = 'mid-PFC';
structRegions(3).electrodes = wholehead_elecs{1,2};
structRegions(4).name = 'right-PFC';
structRegions(4).electrodes = wholehead_elecs{3,2};
structRegions(5).name = 'all-Central';
structRegions(5).electrodes = cat(2, wholehead_elecs{8,2},wholehead_elecs{4,2}, wholehead_elecs{9,2});
structRegions(6).name = 'left-Central';
structRegions(6).electrodes = wholehead_elecs{8,2};
structRegions(7).name = 'mid-Central';
structRegions(7).electrodes = wholehead_elecs{4,2};
structRegions(8).name = 'right-Central';
structRegions(8).electrodes = wholehead_elecs{9,2};
structRegions(9).name = 'all-Post';
structRegions(9).electrodes = cat(2, wholehead_elecs{6,2},wholehead_elecs{5,2}, wholehead_elecs{7,2});
structRegions(10).name = 'left-Post';
structRegions(10).electrodes = wholehead_elecs{6,2};
structRegions(11).name = 'mid-Post';
structRegions(11).electrodes = wholehead_elecs{5,2};
structRegions(12).name = 'right-Post';
structRegions(12).electrodes = wholehead_elecs{7,2};

clear wholehead_elecs

%begin ui algorithm--------------------------------------------------------

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

%prompt for output folder

%create selection list
lstOutput = {'Same as input folder','Custom'};

%show dialog
[idxSelection boolOK] = listdlg('Name','Output directory','PromptString','Select output directory name setting:', ...
                                'ListString',lstOutput,'ListSize',[200 300],'SelectionMode','single');
if boolOK
    switch idxSelection
        case 1  %same as input
            strOutfolder = strInfolder;
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
for i = 1:length(lstInfiles)
    waitbar((i-1)/length(lstInfiles),h,['Loading... ',strrep(lstInfiles{i},'_','\_')]);
    
    load([dirAnalysis,lstInfiles{i}])

    %save into data structure
    dataWavelet(i) = Data;
    if exist('History','var')
        dataHistory(i) = History;
    end
    
    waitbar(i/length(lstInfiles),h);

    clear Data History
end
delete(h);

%cleanup comparisons to only the ones that have required files found
lstDeleteIndexes = [];
for i = 1:length(structComparisons)
    for j = 1:length(structComparisons(i).files)
        if all(cellfun(@isempty,regexp(lstInfiles,structComparisons(i).files(j))))
            %exclude this comparison if any required file is not found
            lstDeleteIndexes(end+1) = i;
        end
    end
end
structComparisons(lstDeleteIndexes) = [];
clear lstDeleteIndexes i j

%loop through all comparisons
figure
for myComparison = structComparisons
    
    %loop through all regions
    for myRegion = structRegions
        %get the electrode list
        lstElectrodes = myRegion.electrodes;
        
        %process both ERSP and ITC
        for cellDatatype = {'ERSP','ITC'}
            strDatatype = cell2mat(cellDatatype);
            
            %find indexes of current data type
            idxDatatype = find(cellfun(@(x) strcmp(x,strDatatype),{dataWavelet.type}));
            
            %check if data is present
            if isempty(idxDatatype);
                %no data is present for current data type, skip
                fprintf('Warning: failed to find %s data for %s %s',strDatatype,myComparison.name,myRegion.name)
                continue
            end
            
            %find the data set from the collection dataWavelet
            %check if manual subtraction is needed (for Qmoney comparisons)
            if (length(myComparison.files) == 4)   %(A-B)-(C-D)
                %parse dosage and condition information from myComparison.files
                strSplit = regexp(myComparison.files{1},'_','split');
                strField1 = strSplit{1};
                strSubfield1 = strSplit{2};
                strSplit = regexp(myComparison.files{2},'_','split');
                strField2 = strSplit{1};
                strSubfield2 = strSplit{2};
                strSplit = regexp(myComparison.files{3},'_','split');
                strField3 = strSplit{1};
                strSubfield3 = strSplit{2};
                strSplit = regexp(myComparison.files{4},'_','split');
                strField4 = strSplit{1};
                strSubfield4 = strSplit{2};
                
                %find where the conditions are stored in the dataWavelet structure
                idxField = find(cellfun(@(x) strcmp(x,strField1),{dataWavelet.field}));
                idxSubfield = find(cellfun(@(x) strcmp(x,strSubfield1),{dataWavelet.subfield}));
                idxGroup1 = intersect(idxField,idxSubfield);
                idxGroup1 = intersect(idxDatatype,idxGroup1);
                idxField = find(cellfun(@(x) strcmp(x,strField2),{dataWavelet.field}));
                idxSubfield = find(cellfun(@(x) strcmp(x,strSubfield2),{dataWavelet.subfield}));
                idxGroup2 = intersect(idxField,idxSubfield);
                idxGroup2 = intersect(idxDatatype,idxGroup2);
                idxField = find(cellfun(@(x) strcmp(x,strField3),{dataWavelet.field}));
                idxSubfield = find(cellfun(@(x) strcmp(x,strSubfield3),{dataWavelet.subfield}));
                idxGroup3 = intersect(idxField,idxSubfield);
                idxGroup3 = intersect(idxDatatype,idxGroup3);
                idxField = find(cellfun(@(x) strcmp(x,strField4),{dataWavelet.field}));
                idxSubfield = find(cellfun(@(x) strcmp(x,strSubfield4),{dataWavelet.subfield}));
                idxGroup4 = intersect(idxField,idxSubfield);
                idxGroup4 = intersect(idxDatatype,idxGroup4);
                
                %perform subtraction, and modify structure of Data
                myData = dataWavelet(idxGroup1);
                myData.field = myComparison.name;
                myData.subfield = '';   %this is ignored
                myData.matrix = (dataWavelet(idxGroup1).matrix - dataWavelet(idxGroup2).matrix) - ...
                                (dataWavelet(idxGroup3).matrix - dataWavelet(idxGroup4).matrix);
                            
                %build source string
                if exist('dataHistory','var')
                    strSource1_name = dataHistory(idxGroup1).source_name;
                    strSource1_id = dataHistory(idxGroup1).source_id;
                    strSource2_name = dataHistory(idxGroup2).source_name;
                    strSource2_id = dataHistory(idxGroup2).source_id;
                    strSource3_name = dataHistory(idxGroup3).source_name;
                    strSource3_id = dataHistory(idxGroup3).source_id;
                    strSource4_name = dataHistory(idxGroup4).source_name;
                    strSource4_id = dataHistory(idxGroup4).source_id;
                    if (strcmp(strSource1_name,strSource2_name) && strcmp(strSource1_id,strSource2_id) && ...
                        strcmp(strSource1_name,strSource3_name) && strcmp(strSource1_id,strSource3_id) && ...
                        strcmp(strSource1_name,strSource4_name) && strcmp(strSource1_id,strSource4_id))
                        strSource = ['Source: ',strrep(strSource1_name,'_','\_'),' (',strSource1_id,')'];
                    elseif (strcmp(strSource1_name,strSource2_name) && strcmp(strSource1_id,strSource2_id) && ...
                            strcmp(strSource3_name,strSource4_name) && strcmp(strSource3_id,strSource4_id))
                        strSource = {['Source: ',strrep(strSource1_name,'_','\_'),' (',strSource1_id,')'],[strrep(strSource3_name,'_','\_'),' (',strSource3_id,')']};
                    else
                        strSource = {['Source: (',strSource1_id,') (',strSource2_id,')'],['(',strSource3_id,') (',strSource4_id,')']};
                    end
                else
                    strSource = ['Source: ',strrep(myData.source,'_','\_')];
                end
                
                clear strField1 strField2 strField3 strField4 strSubfield1 strSubfield2 strSubfield3 strSubfield4
                clear idxField idxSubfield idxGroup1 idxGroup2 idxGroup3 idxGroup4
                clear strSource1_name strSource1_id strSource2_name strSource2_id strSource3_name strSource3_id strSource4_name strSource4_id
            elseif (length(myComparison.files) == 2)   %(A-B)
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
            elseif (length(myComparison.files) == 1)   %(A-Null)
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
            else
                error('Invalid number of datasets found (should be 1,2,4): %i',length(myComparison.files))
            end
                
            %check if subject exclusion is enabled
            if boolUseSubjectExclusion
                idxDelete = find(ismember(myData.subject_list,lstExcludedSubjects));
                myData.matrix(idxDelete,:,:,:) = [];
                
                clear idxDelete
            end
            
            %extract electrode data
            myData.matrix = myData.matrix(:,lstElectrodes,:,:);
            
            %average down the region to 1 giant electrode <subject X frequency X time>
            myData.matrix = squeeze(nanmean(myData.matrix,2));
            
            %average down all subjects to get group level results <frequency X time>
            myData.matrix = squeeze(nanmean(myData.matrix,1));
            
            %set color scale
            if strcmp(strDatatype,'ERSP')
                idxColorScale = find(cellfun(@(x) strcmp(x,'ERSP'),{structColorScale.type}));
            elseif (min(myData.matrix(:)) >= 0)
                idxColorScale = find(cellfun(@(x) strcmp(x,'ITC single condition'),{structColorScale.type}));
            else
                idxColorScale = find(cellfun(@(x) strcmp(x,'ITC contrasts'),{structColorScale.type}));
            end
            
            %build filename
            strBuilder = [strFilePrefix,'_'];
            strBuilder = [strBuilder,strDatatype,'_'];
            strBuilder = [strBuilder,myComparison.name,'_'];
            strBuilder = [strBuilder,myRegion.name];
            
            %build title string
            if boolPublicationMode
                strTitle = ' ';
            else
                strTitle = ['{\bf{\color{blue}',strDatatype,'}} ',myComparison.name,' {\color{red}',myRegion.name,'}'];
            end
            
            %graph the result
            clf
            h = tftopo(myData.matrix,myData.time,myData.frequency,'title',strTitle,'verbose','off');
            if boolUseCustomScaling
                caxis(structColorScale(idxColorScale).limits), ...
                colorbar;
            else
                colorbar;
            end
            if ~boolPublicationMode
                lblText = axes('Units','Normal','Position',[0.35 -0.86 .85 .85],'Visible','off');
                set(get(lblText,'Title'),'Visible','on'), title(strSource,'FontSize',6);
            end
            saveas(gcf,[dirOutput,strBuilder,'.tif'],'tiff');
            
            clear myData idxDatatype strDatatype lblText idxColorScale strBuilder strTitle strSource
        end
        clear lstElectrodes
    end
end

close,beep