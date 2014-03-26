%Extract statistics from the behave performance data
%Author: Julian Y. Cheng
%
%Note: this loads the CSV files converted from the EDAT text files, not the
%      MAT files because apparantly all non-hit trials are absent
%
%Changelog:
%   3/13/2013: re-converted text files to include block numbers. The
%              indexing used in the algorithm now reflect the new file
%              format. Added per-block processing to algorithm. Also,
%              modified to report excluded trial count on a per-block
%              basis.

clc,clear

%user-defined vars
flagDebug = 0;  %avoids writing to file

%task-specific vars
intMinRT = 200; %any value outside of [min,max] will be excluded, units in ms
intMaxRT = 1200;
dirRoot = '/nfs/erp-modaf/mc/MC-epoched-Feb11-Glenn/RT_analysis/';
dirDrug = [dirRoot,'Drug/'];
dirPlac = [dirRoot,'Plac/'];
dirOutput = dirRoot;
strOutfile = 'MCPC_behave_stats.csv';
strGreenLabel = 'green';
strRedLabel = 'red';
arrHeader = ... %output header; Note: if you change this, you need to change how the data is stored
    [{'Subject'} ...
     {'PGQ_ACC_b1'} {'PGQ_ACC_b2'} {'PGQ_ACC_b3'} {'PGQ_ACC_b4'} {'PGQ_ACC_b5'} {'PGQ_ACC_b6'} {'PGQ_ACC_b7'} {'PGQ_ACC_b8'} {'PGQ_ACC_all'} ...
     {'PRQ_ACC_b1'} {'PRQ_ACC_b2'} {'PRQ_ACC_b3'} {'PRQ_ACC_b4'} {'PRQ_ACC_b5'} {'PRQ_ACC_b6'} {'PRQ_ACC_b7'} {'PRQ_ACC_b8'} {'PRQ_ACC_all'} ...
     {'DGQ_ACC_b1'} {'DGQ_ACC_b2'} {'DGQ_ACC_b3'} {'DGQ_ACC_b4'} {'DGQ_ACC_b5'} {'DGQ_ACC_b6'} {'DGQ_ACC_b7'} {'DGQ_ACC_b8'} {'DGQ_ACC_all'} ...
     {'DRQ_ACC_b1'} {'DRQ_ACC_b2'} {'DRQ_ACC_b3'} {'DRQ_ACC_b4'} {'DRQ_ACC_b5'} {'DRQ_ACC_b6'} {'DRQ_ACC_b7'} {'DRQ_ACC_b8'} {'DRQ_ACC_all'} ...
     {'PGQ_RT_b1'} {'PGQ_RT_b2'} {'PGQ_RT_b3'} {'PGQ_RT_b4'} {'PGQ_RT_b5'} {'PGQ_RT_b6'} {'PGQ_RT_b7'} {'PGQ_RT_b8'} {'PGQ_RT_all'} ...
     {'PRQ_RT_b1'} {'PRQ_RT_b2'} {'PRQ_RT_b3'} {'PRQ_RT_b4'} {'PRQ_RT_b5'} {'PRQ_RT_b6'} {'PRQ_RT_b7'} {'PRQ_RT_b8'} {'PRQ_RT_all'} ...
     {'DGQ_RT_b1'} {'DGQ_RT_b2'} {'DGQ_RT_b3'} {'DGQ_RT_b4'} {'DGQ_RT_b5'} {'DGQ_RT_b6'} {'DGQ_RT_b7'} {'DGQ_RT_b8'} {'DGQ_RT_all'} ...
     {'DRQ_RT_b1'} {'DRQ_RT_b2'} {'DRQ_RT_b3'} {'DRQ_RT_b4'} {'DRQ_RT_b5'} {'DRQ_RT_b6'} {'DRQ_RT_b7'} {'DRQ_RT_b8'} {'DRQ_RT_all'} ...
     {'PGQ_RT-SD_b1'} {'PGQ_RT-SD_b2'} {'PGQ_RT-SD_b3'} {'PGQ_RT-SD_b4'} {'PGQ_RT-SD_b5'} {'PGQ_RT-SD_b6'} {'PGQ_RT-SD_b7'} {'PGQ_RT-SD_b8'} {'PGQ_RT-SD_all'} ...
     {'PRQ_RT-SD_b1'} {'PRQ_RT-SD_b2'} {'PRQ_RT-SD_b3'} {'PRQ_RT-SD_b4'} {'PRQ_RT-SD_b5'} {'PRQ_RT-SD_b6'} {'PRQ_RT-SD_b7'} {'PRQ_RT-SD_b8'} {'PRQ_RT-SD_all'} ...
     {'DGQ_RT-SD_b1'} {'DGQ_RT-SD_b2'} {'DGQ_RT-SD_b3'} {'DGQ_RT-SD_b4'} {'DGQ_RT-SD_b5'} {'DGQ_RT-SD_b6'} {'DGQ_RT-SD_b7'} {'DGQ_RT-SD_b8'} {'DGQ_RT-SD_all'} ...
     {'DRQ_RT-SD_b1'} {'DRQ_RT-SD_b2'} {'DRQ_RT-SD_b3'} {'DRQ_RT-SD_b4'} {'DRQ_RT-SD_b5'} {'DRQ_RT-SD_b6'} {'DRQ_RT-SD_b7'} {'DRQ_RT-SD_b8'} {'DRQ_RT-SD_all'} ...
     {'EXC_Plac'} {'%EXC_P'} {'EXC_Drug'} {'%EXC_D'}];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%begin algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Behave performance extraction for statistical analysis\nAuthor: Julian Y. Cheng\n\n');

%get list of drug and plac files
lstDrugInfiles = dir([dirDrug,'*.csv']);
lstDrugInfiles = {lstDrugInfiles.name};
lstPlacInfiles = dir([dirPlac,'*.csv']);
lstPlacInfiles = {lstPlacInfiles.name};

%check if there are files to process
if isempty(lstDrugInfiles)
    error('Files not found for drug: %s',dirDrug)
elseif isempty(lstPlacInfiles)
    error('Files not found for plac: %s',dirPlac)
end

%process both drug and plac data files
%loads in 1 at a time and adds to summary data matrix

matrixAllSubjectData = NaN(max(length(lstDrugInfiles),length(lstPlacInfiles)),length(arrHeader));

fprintf('Processing drug MAT files:\n')
for cellInfile = lstDrugInfiles     %all drug files
    strInfile = cell2mat(cellInfile);
    
    %get the subject id from the filename
    lstStringSplit = regexp(strInfile,'[.]','split');
    strSubjectID = lstStringSplit{1};
    
    fprintf('Subject %s: \n',strSubjectID)

    %load data and condense to 1 variable
    structData = importdata([dirDrug,strInfile]);
    behave = structData.textdata;
    behave(2:end,3:end) = num2cell(structData.data);
    behave(1,:) = [];   %delete header
    behave(:,1) = num2cell(str2num(cell2mat(behave(:,1)))); %convert string block numbers into integer block numbers
    
    %find all entries for each block
    lstBlockIndexes1 = find(cellfun(@(x) isequal(x,1),behave(:,1)));
    lstBlockIndexes2 = find(cellfun(@(x) isequal(x,2),behave(:,1)));
    lstBlockIndexes3 = find(cellfun(@(x) isequal(x,3),behave(:,1)));
    lstBlockIndexes4 = find(cellfun(@(x) isequal(x,4),behave(:,1)));
    lstBlockIndexes5 = find(cellfun(@(x) isequal(x,5),behave(:,1)));
    lstBlockIndexes6 = find(cellfun(@(x) isequal(x,6),behave(:,1)));
    lstBlockIndexes7 = find(cellfun(@(x) isequal(x,7),behave(:,1)));
    lstBlockIndexes8 = find(cellfun(@(x) isequal(x,8),behave(:,1)));
    
    %check exclusion criteria
    lstTooFastIndexes = find(cell2mat(behave(:,4)) < intMinRT);
    lstTooSlowIndexes = find(cell2mat(behave(:,4)) > intMaxRT);
    fprintf('\tBlock 1: %2i (too fast) %2i (too slow)\n',length(intersect(lstTooFastIndexes,lstBlockIndexes1)),length(intersect(lstTooSlowIndexes,lstBlockIndexes1)))
    fprintf('\tBlock 2: %2i (too fast) %2i (too slow)\n',length(intersect(lstTooFastIndexes,lstBlockIndexes2)),length(intersect(lstTooSlowIndexes,lstBlockIndexes2)))
    fprintf('\tBlock 3: %2i (too fast) %2i (too slow)\n',length(intersect(lstTooFastIndexes,lstBlockIndexes3)),length(intersect(lstTooSlowIndexes,lstBlockIndexes3)))
    fprintf('\tBlock 4: %2i (too fast) %2i (too slow)\n',length(intersect(lstTooFastIndexes,lstBlockIndexes4)),length(intersect(lstTooSlowIndexes,lstBlockIndexes4)))
    fprintf('\tBlock 5: %2i (too fast) %2i (too slow)\n',length(intersect(lstTooFastIndexes,lstBlockIndexes5)),length(intersect(lstTooSlowIndexes,lstBlockIndexes5)))
    fprintf('\tBlock 6: %2i (too fast) %2i (too slow)\n',length(intersect(lstTooFastIndexes,lstBlockIndexes6)),length(intersect(lstTooSlowIndexes,lstBlockIndexes6)))
    fprintf('\tBlock 7: %2i (too fast) %2i (too slow)\n',length(intersect(lstTooFastIndexes,lstBlockIndexes7)),length(intersect(lstTooSlowIndexes,lstBlockIndexes7)))
    fprintf('\tBlock 8: %2i (too fast) %2i (too slow)\n',length(intersect(lstTooFastIndexes,lstBlockIndexes8)),length(intersect(lstTooSlowIndexes,lstBlockIndexes8)))
    lstExcludedIndexes = [lstTooFastIndexes;lstTooSlowIndexes];
    lstAcceptedIndexes = setdiff(1:size(behave,1),lstExcludedIndexes);
    
    %check if block lengths are the same
    if (length(unique([length(lstBlockIndexes1),length(lstBlockIndexes2),length(lstBlockIndexes3),length(lstBlockIndexes4),length(lstBlockIndexes5),length(lstBlockIndexes6),length(lstBlockIndexes7),length(lstBlockIndexes8)])) ~= 1)
        fprintf('Warning: block lengths are NOT the same\n')
    end
    
    %find all entries for each condition
    lstGreenIndexes = find(cellfun(@(x) strcmp(x,strGreenLabel),behave(lstAcceptedIndexes,2)));
    lstRedIndexes = find(cellfun(@(x) strcmp(x,strRedLabel),behave(lstAcceptedIndexes,2)));
    
    %get statistical measures
    dblGreenACC_all = sum(cell2mat(behave(lstGreenIndexes,3))) / length(lstGreenIndexes);
    dblGreenACC_blk1 = sum(cell2mat(behave(intersect(lstGreenIndexes,lstBlockIndexes1),3))) / length(intersect(lstGreenIndexes,lstBlockIndexes1));
    dblGreenACC_blk2 = sum(cell2mat(behave(intersect(lstGreenIndexes,lstBlockIndexes2),3))) / length(intersect(lstGreenIndexes,lstBlockIndexes2));
    dblGreenACC_blk3 = sum(cell2mat(behave(intersect(lstGreenIndexes,lstBlockIndexes3),3))) / length(intersect(lstGreenIndexes,lstBlockIndexes3));
    dblGreenACC_blk4 = sum(cell2mat(behave(intersect(lstGreenIndexes,lstBlockIndexes4),3))) / length(intersect(lstGreenIndexes,lstBlockIndexes4));
    dblGreenACC_blk5 = sum(cell2mat(behave(intersect(lstGreenIndexes,lstBlockIndexes5),3))) / length(intersect(lstGreenIndexes,lstBlockIndexes5));
    dblGreenACC_blk6 = sum(cell2mat(behave(intersect(lstGreenIndexes,lstBlockIndexes6),3))) / length(intersect(lstGreenIndexes,lstBlockIndexes6));
    dblGreenACC_blk7 = sum(cell2mat(behave(intersect(lstGreenIndexes,lstBlockIndexes7),3))) / length(intersect(lstGreenIndexes,lstBlockIndexes7));
    dblGreenACC_blk8 = sum(cell2mat(behave(intersect(lstGreenIndexes,lstBlockIndexes8),3))) / length(intersect(lstGreenIndexes,lstBlockIndexes8));
    dblRedACC_all = sum(cell2mat(behave(lstRedIndexes,3))) / length(lstRedIndexes);
    dblRedACC_blk1 = sum(cell2mat(behave(intersect(lstRedIndexes,lstBlockIndexes1),3))) / length(intersect(lstRedIndexes,lstBlockIndexes1));
    dblRedACC_blk2 = sum(cell2mat(behave(intersect(lstRedIndexes,lstBlockIndexes2),3))) / length(intersect(lstRedIndexes,lstBlockIndexes2));
    dblRedACC_blk3 = sum(cell2mat(behave(intersect(lstRedIndexes,lstBlockIndexes3),3))) / length(intersect(lstRedIndexes,lstBlockIndexes3));
    dblRedACC_blk4 = sum(cell2mat(behave(intersect(lstRedIndexes,lstBlockIndexes4),3))) / length(intersect(lstRedIndexes,lstBlockIndexes4));
    dblRedACC_blk5 = sum(cell2mat(behave(intersect(lstRedIndexes,lstBlockIndexes5),3))) / length(intersect(lstRedIndexes,lstBlockIndexes5));
    dblRedACC_blk6 = sum(cell2mat(behave(intersect(lstRedIndexes,lstBlockIndexes6),3))) / length(intersect(lstRedIndexes,lstBlockIndexes6));
    dblRedACC_blk7 = sum(cell2mat(behave(intersect(lstRedIndexes,lstBlockIndexes7),3))) / length(intersect(lstRedIndexes,lstBlockIndexes7));
    dblRedACC_blk8 = sum(cell2mat(behave(intersect(lstRedIndexes,lstBlockIndexes8),3))) / length(intersect(lstRedIndexes,lstBlockIndexes8));
    dblGreenRT_all = mean(cell2mat(behave(lstGreenIndexes,4)));
    dblGreenRT_blk1 = mean(cell2mat(behave(intersect(lstGreenIndexes,lstBlockIndexes1),4)));
    dblGreenRT_blk2 = mean(cell2mat(behave(intersect(lstGreenIndexes,lstBlockIndexes2),4)));
    dblGreenRT_blk3 = mean(cell2mat(behave(intersect(lstGreenIndexes,lstBlockIndexes3),4)));
    dblGreenRT_blk4 = mean(cell2mat(behave(intersect(lstGreenIndexes,lstBlockIndexes4),4)));
    dblGreenRT_blk5 = mean(cell2mat(behave(intersect(lstGreenIndexes,lstBlockIndexes5),4)));
    dblGreenRT_blk6 = mean(cell2mat(behave(intersect(lstGreenIndexes,lstBlockIndexes6),4)));
    dblGreenRT_blk7 = mean(cell2mat(behave(intersect(lstGreenIndexes,lstBlockIndexes7),4)));
    dblGreenRT_blk8 = mean(cell2mat(behave(intersect(lstGreenIndexes,lstBlockIndexes8),4)));
    dblRedRT_all = mean(cell2mat(behave(lstRedIndexes,4)));
    dblRedRT_blk1 = mean(cell2mat(behave(intersect(lstRedIndexes,lstBlockIndexes1),4)));
    dblRedRT_blk2 = mean(cell2mat(behave(intersect(lstRedIndexes,lstBlockIndexes2),4)));
    dblRedRT_blk3 = mean(cell2mat(behave(intersect(lstRedIndexes,lstBlockIndexes3),4)));
    dblRedRT_blk4 = mean(cell2mat(behave(intersect(lstRedIndexes,lstBlockIndexes4),4)));
    dblRedRT_blk5 = mean(cell2mat(behave(intersect(lstRedIndexes,lstBlockIndexes5),4)));
    dblRedRT_blk6 = mean(cell2mat(behave(intersect(lstRedIndexes,lstBlockIndexes6),4)));
    dblRedRT_blk7 = mean(cell2mat(behave(intersect(lstRedIndexes,lstBlockIndexes7),4)));
    dblRedRT_blk8 = mean(cell2mat(behave(intersect(lstRedIndexes,lstBlockIndexes8),4)));
    dblGreenRTSD_all = std(cell2mat(behave(lstGreenIndexes,4)));
    dblGreenRTSD_blk1 = std(cell2mat(behave(intersect(lstGreenIndexes,lstBlockIndexes1),4)));
    dblGreenRTSD_blk2 = std(cell2mat(behave(intersect(lstGreenIndexes,lstBlockIndexes2),4)));
    dblGreenRTSD_blk3 = std(cell2mat(behave(intersect(lstGreenIndexes,lstBlockIndexes3),4)));
    dblGreenRTSD_blk4 = std(cell2mat(behave(intersect(lstGreenIndexes,lstBlockIndexes4),4)));
    dblGreenRTSD_blk5 = std(cell2mat(behave(intersect(lstGreenIndexes,lstBlockIndexes5),4)));
    dblGreenRTSD_blk6 = std(cell2mat(behave(intersect(lstGreenIndexes,lstBlockIndexes6),4)));
    dblGreenRTSD_blk7 = std(cell2mat(behave(intersect(lstGreenIndexes,lstBlockIndexes7),4)));
    dblGreenRTSD_blk8 = std(cell2mat(behave(intersect(lstGreenIndexes,lstBlockIndexes8),4)));
    dblRedRTSD_all = std(cell2mat(behave(lstRedIndexes,4)));
    dblRedRTSD_blk1 = std(cell2mat(behave(intersect(lstRedIndexes,lstBlockIndexes1),4)));
    dblRedRTSD_blk2 = std(cell2mat(behave(intersect(lstRedIndexes,lstBlockIndexes2),4)));
    dblRedRTSD_blk3 = std(cell2mat(behave(intersect(lstRedIndexes,lstBlockIndexes3),4)));
    dblRedRTSD_blk4 = std(cell2mat(behave(intersect(lstRedIndexes,lstBlockIndexes4),4)));
    dblRedRTSD_blk5 = std(cell2mat(behave(intersect(lstRedIndexes,lstBlockIndexes5),4)));
    dblRedRTSD_blk6 = std(cell2mat(behave(intersect(lstRedIndexes,lstBlockIndexes6),4)));
    dblRedRTSD_blk7 = std(cell2mat(behave(intersect(lstRedIndexes,lstBlockIndexes7),4)));
    dblRedRTSD_blk8 = std(cell2mat(behave(intersect(lstRedIndexes,lstBlockIndexes8),4)));
    
    %save into all subject data matrix
    idxNextEntry = find(isnan(matrixAllSubjectData(:,1)),1);
    matrixAllSubjectData(idxNextEntry,1) = str2double(strSubjectID);
    matrixAllSubjectData(idxNextEntry,20) = dblGreenACC_blk1;
    matrixAllSubjectData(idxNextEntry,21) = dblGreenACC_blk2;
    matrixAllSubjectData(idxNextEntry,22) = dblGreenACC_blk3;
    matrixAllSubjectData(idxNextEntry,23) = dblGreenACC_blk4;
    matrixAllSubjectData(idxNextEntry,24) = dblGreenACC_blk5;
    matrixAllSubjectData(idxNextEntry,25) = dblGreenACC_blk6;
    matrixAllSubjectData(idxNextEntry,26) = dblGreenACC_blk7;
    matrixAllSubjectData(idxNextEntry,27) = dblGreenACC_blk8;
    matrixAllSubjectData(idxNextEntry,28) = dblGreenACC_all;
    matrixAllSubjectData(idxNextEntry,29) = dblRedACC_blk1;
    matrixAllSubjectData(idxNextEntry,30) = dblRedACC_blk2;
    matrixAllSubjectData(idxNextEntry,31) = dblRedACC_blk3;
    matrixAllSubjectData(idxNextEntry,32) = dblRedACC_blk4;
    matrixAllSubjectData(idxNextEntry,33) = dblRedACC_blk5;
    matrixAllSubjectData(idxNextEntry,34) = dblRedACC_blk6;
    matrixAllSubjectData(idxNextEntry,35) = dblRedACC_blk7;
    matrixAllSubjectData(idxNextEntry,36) = dblRedACC_blk8;
    matrixAllSubjectData(idxNextEntry,37) = dblRedACC_all;
    matrixAllSubjectData(idxNextEntry,56) = dblGreenRT_blk1;
    matrixAllSubjectData(idxNextEntry,57) = dblGreenRT_blk2;
    matrixAllSubjectData(idxNextEntry,58) = dblGreenRT_blk3;
    matrixAllSubjectData(idxNextEntry,59) = dblGreenRT_blk4;
    matrixAllSubjectData(idxNextEntry,60) = dblGreenRT_blk5;
    matrixAllSubjectData(idxNextEntry,61) = dblGreenRT_blk6;
    matrixAllSubjectData(idxNextEntry,62) = dblGreenRT_blk7;
    matrixAllSubjectData(idxNextEntry,63) = dblGreenRT_blk8;
    matrixAllSubjectData(idxNextEntry,64) = dblGreenRT_all;
    matrixAllSubjectData(idxNextEntry,65) = dblRedRT_blk1;
    matrixAllSubjectData(idxNextEntry,66) = dblRedRT_blk2;
    matrixAllSubjectData(idxNextEntry,67) = dblRedRT_blk3;
    matrixAllSubjectData(idxNextEntry,68) = dblRedRT_blk4;
    matrixAllSubjectData(idxNextEntry,69) = dblRedRT_blk5;
    matrixAllSubjectData(idxNextEntry,70) = dblRedRT_blk6;
    matrixAllSubjectData(idxNextEntry,71) = dblRedRT_blk7;
    matrixAllSubjectData(idxNextEntry,72) = dblRedRT_blk8;
    matrixAllSubjectData(idxNextEntry,73) = dblRedRT_all;
    matrixAllSubjectData(idxNextEntry,92) = dblGreenRTSD_blk1;
    matrixAllSubjectData(idxNextEntry,93) = dblGreenRTSD_blk2;
    matrixAllSubjectData(idxNextEntry,94) = dblGreenRTSD_blk3;
    matrixAllSubjectData(idxNextEntry,95) = dblGreenRTSD_blk4;
    matrixAllSubjectData(idxNextEntry,96) = dblGreenRTSD_blk5;
    matrixAllSubjectData(idxNextEntry,97) = dblGreenRTSD_blk6;
    matrixAllSubjectData(idxNextEntry,98) = dblGreenRTSD_blk7;
    matrixAllSubjectData(idxNextEntry,99) = dblGreenRTSD_blk8;
    matrixAllSubjectData(idxNextEntry,100) = dblGreenRTSD_all;
    matrixAllSubjectData(idxNextEntry,101) = dblRedRTSD_blk1;
    matrixAllSubjectData(idxNextEntry,102) = dblRedRTSD_blk2;
    matrixAllSubjectData(idxNextEntry,103) = dblRedRTSD_blk3;
    matrixAllSubjectData(idxNextEntry,104) = dblRedRTSD_blk4;
    matrixAllSubjectData(idxNextEntry,105) = dblRedRTSD_blk5;
    matrixAllSubjectData(idxNextEntry,106) = dblRedRTSD_blk6;
    matrixAllSubjectData(idxNextEntry,107) = dblRedRTSD_blk7;
    matrixAllSubjectData(idxNextEntry,108) = dblRedRTSD_blk8;
    matrixAllSubjectData(idxNextEntry,109) = dblRedRTSD_all;
    matrixAllSubjectData(idxNextEntry,112) = length(lstExcludedIndexes);
    matrixAllSubjectData(idxNextEntry,113) = length(lstExcludedIndexes)/size(behave,1)*100;
    
    clear behave structData
end
fprintf('Processing plac MAT files:\n')
for cellInfile = lstPlacInfiles     %all plac files
    strInfile = cell2mat(cellInfile);
    
    %get the subject id from the filename
    lstStringSplit = regexp(strInfile,'[.]','split');
    strSubjectID = lstStringSplit{1};
    
    fprintf('Subject %s:\n',strSubjectID)
    
    %load data and condense to 1 variable
    structData = importdata([dirPlac,strInfile]);
    behave = structData.textdata;
    behave(2:end,3:end) = num2cell(structData.data);
    behave(1,:) = [];   %delete header
    behave(:,1) = num2cell(str2num(cell2mat(behave(:,1)))); %convert string block numbers into integer block numbers
    
    %find all entries for each block
    lstBlockIndexes1 = find(cellfun(@(x) isequal(x,1),behave(:,1)));
    lstBlockIndexes2 = find(cellfun(@(x) isequal(x,2),behave(:,1)));
    lstBlockIndexes3 = find(cellfun(@(x) isequal(x,3),behave(:,1)));
    lstBlockIndexes4 = find(cellfun(@(x) isequal(x,4),behave(:,1)));
    lstBlockIndexes5 = find(cellfun(@(x) isequal(x,5),behave(:,1)));
    lstBlockIndexes6 = find(cellfun(@(x) isequal(x,6),behave(:,1)));
    lstBlockIndexes7 = find(cellfun(@(x) isequal(x,7),behave(:,1)));
    lstBlockIndexes8 = find(cellfun(@(x) isequal(x,8),behave(:,1)));
    
    %check exclusion criteria
    lstTooFastIndexes = find(cell2mat(behave(:,4)) < intMinRT);
    lstTooSlowIndexes = find(cell2mat(behave(:,4)) > intMaxRT);
    fprintf('\tBlock 1: %2i (too fast) %2i (too slow)\n',length(intersect(lstTooFastIndexes,lstBlockIndexes1)),length(intersect(lstTooSlowIndexes,lstBlockIndexes1)))
    fprintf('\tBlock 2: %2i (too fast) %2i (too slow)\n',length(intersect(lstTooFastIndexes,lstBlockIndexes2)),length(intersect(lstTooSlowIndexes,lstBlockIndexes2)))
    fprintf('\tBlock 3: %2i (too fast) %2i (too slow)\n',length(intersect(lstTooFastIndexes,lstBlockIndexes3)),length(intersect(lstTooSlowIndexes,lstBlockIndexes3)))
    fprintf('\tBlock 4: %2i (too fast) %2i (too slow)\n',length(intersect(lstTooFastIndexes,lstBlockIndexes4)),length(intersect(lstTooSlowIndexes,lstBlockIndexes4)))
    fprintf('\tBlock 5: %2i (too fast) %2i (too slow)\n',length(intersect(lstTooFastIndexes,lstBlockIndexes5)),length(intersect(lstTooSlowIndexes,lstBlockIndexes5)))
    fprintf('\tBlock 6: %2i (too fast) %2i (too slow)\n',length(intersect(lstTooFastIndexes,lstBlockIndexes6)),length(intersect(lstTooSlowIndexes,lstBlockIndexes6)))
    fprintf('\tBlock 7: %2i (too fast) %2i (too slow)\n',length(intersect(lstTooFastIndexes,lstBlockIndexes7)),length(intersect(lstTooSlowIndexes,lstBlockIndexes7)))
    fprintf('\tBlock 8: %2i (too fast) %2i (too slow)\n',length(intersect(lstTooFastIndexes,lstBlockIndexes8)),length(intersect(lstTooSlowIndexes,lstBlockIndexes8)))
    lstExcludedIndexes = [lstTooFastIndexes;lstTooSlowIndexes];
    lstAcceptedIndexes = setdiff(1:size(behave,1),lstExcludedIndexes);
    
    %check if block lengths are the same
    if (length(unique([length(lstBlockIndexes1),length(lstBlockIndexes2),length(lstBlockIndexes3),length(lstBlockIndexes4),length(lstBlockIndexes5),length(lstBlockIndexes6),length(lstBlockIndexes7),length(lstBlockIndexes8)])) ~= 1)
        fprintf('Warning: block lengths are NOT the same\n')
    end
    
    %find all entries for each condition
    lstGreenIndexes = find(cellfun(@(x) strcmp(x,strGreenLabel),behave(lstAcceptedIndexes,2)));
    lstRedIndexes = find(cellfun(@(x) strcmp(x,strRedLabel),behave(lstAcceptedIndexes,2)));
    
    %get statistical measures
    dblGreenACC_all = sum(cell2mat(behave(lstGreenIndexes,3))) / length(lstGreenIndexes);
    dblGreenACC_blk1 = sum(cell2mat(behave(intersect(lstGreenIndexes,lstBlockIndexes1),3))) / length(intersect(lstGreenIndexes,lstBlockIndexes1));
    dblGreenACC_blk2 = sum(cell2mat(behave(intersect(lstGreenIndexes,lstBlockIndexes2),3))) / length(intersect(lstGreenIndexes,lstBlockIndexes2));
    dblGreenACC_blk3 = sum(cell2mat(behave(intersect(lstGreenIndexes,lstBlockIndexes3),3))) / length(intersect(lstGreenIndexes,lstBlockIndexes3));
    dblGreenACC_blk4 = sum(cell2mat(behave(intersect(lstGreenIndexes,lstBlockIndexes4),3))) / length(intersect(lstGreenIndexes,lstBlockIndexes4));
    dblGreenACC_blk5 = sum(cell2mat(behave(intersect(lstGreenIndexes,lstBlockIndexes5),3))) / length(intersect(lstGreenIndexes,lstBlockIndexes5));
    dblGreenACC_blk6 = sum(cell2mat(behave(intersect(lstGreenIndexes,lstBlockIndexes6),3))) / length(intersect(lstGreenIndexes,lstBlockIndexes6));
    dblGreenACC_blk7 = sum(cell2mat(behave(intersect(lstGreenIndexes,lstBlockIndexes7),3))) / length(intersect(lstGreenIndexes,lstBlockIndexes7));
    dblGreenACC_blk8 = sum(cell2mat(behave(intersect(lstGreenIndexes,lstBlockIndexes8),3))) / length(intersect(lstGreenIndexes,lstBlockIndexes8));
    dblRedACC_all = sum(cell2mat(behave(lstRedIndexes,3))) / length(lstRedIndexes);
    dblRedACC_blk1 = sum(cell2mat(behave(intersect(lstRedIndexes,lstBlockIndexes1),3))) / length(intersect(lstRedIndexes,lstBlockIndexes1));
    dblRedACC_blk2 = sum(cell2mat(behave(intersect(lstRedIndexes,lstBlockIndexes2),3))) / length(intersect(lstRedIndexes,lstBlockIndexes2));
    dblRedACC_blk3 = sum(cell2mat(behave(intersect(lstRedIndexes,lstBlockIndexes3),3))) / length(intersect(lstRedIndexes,lstBlockIndexes3));
    dblRedACC_blk4 = sum(cell2mat(behave(intersect(lstRedIndexes,lstBlockIndexes4),3))) / length(intersect(lstRedIndexes,lstBlockIndexes4));
    dblRedACC_blk5 = sum(cell2mat(behave(intersect(lstRedIndexes,lstBlockIndexes5),3))) / length(intersect(lstRedIndexes,lstBlockIndexes5));
    dblRedACC_blk6 = sum(cell2mat(behave(intersect(lstRedIndexes,lstBlockIndexes6),3))) / length(intersect(lstRedIndexes,lstBlockIndexes6));
    dblRedACC_blk7 = sum(cell2mat(behave(intersect(lstRedIndexes,lstBlockIndexes7),3))) / length(intersect(lstRedIndexes,lstBlockIndexes7));
    dblRedACC_blk8 = sum(cell2mat(behave(intersect(lstRedIndexes,lstBlockIndexes8),3))) / length(intersect(lstRedIndexes,lstBlockIndexes8));
    dblGreenRT_all = mean(cell2mat(behave(lstGreenIndexes,4)));
    dblGreenRT_blk1 = mean(cell2mat(behave(intersect(lstGreenIndexes,lstBlockIndexes1),4)));
    dblGreenRT_blk2 = mean(cell2mat(behave(intersect(lstGreenIndexes,lstBlockIndexes2),4)));
    dblGreenRT_blk3 = mean(cell2mat(behave(intersect(lstGreenIndexes,lstBlockIndexes3),4)));
    dblGreenRT_blk4 = mean(cell2mat(behave(intersect(lstGreenIndexes,lstBlockIndexes4),4)));
    dblGreenRT_blk5 = mean(cell2mat(behave(intersect(lstGreenIndexes,lstBlockIndexes5),4)));
    dblGreenRT_blk6 = mean(cell2mat(behave(intersect(lstGreenIndexes,lstBlockIndexes6),4)));
    dblGreenRT_blk7 = mean(cell2mat(behave(intersect(lstGreenIndexes,lstBlockIndexes7),4)));
    dblGreenRT_blk8 = mean(cell2mat(behave(intersect(lstGreenIndexes,lstBlockIndexes8),4)));
    dblRedRT_all = mean(cell2mat(behave(lstRedIndexes,4)));
    dblRedRT_blk1 = mean(cell2mat(behave(intersect(lstRedIndexes,lstBlockIndexes1),4)));
    dblRedRT_blk2 = mean(cell2mat(behave(intersect(lstRedIndexes,lstBlockIndexes2),4)));
    dblRedRT_blk3 = mean(cell2mat(behave(intersect(lstRedIndexes,lstBlockIndexes3),4)));
    dblRedRT_blk4 = mean(cell2mat(behave(intersect(lstRedIndexes,lstBlockIndexes4),4)));
    dblRedRT_blk5 = mean(cell2mat(behave(intersect(lstRedIndexes,lstBlockIndexes5),4)));
    dblRedRT_blk6 = mean(cell2mat(behave(intersect(lstRedIndexes,lstBlockIndexes6),4)));
    dblRedRT_blk7 = mean(cell2mat(behave(intersect(lstRedIndexes,lstBlockIndexes7),4)));
    dblRedRT_blk8 = mean(cell2mat(behave(intersect(lstRedIndexes,lstBlockIndexes8),4)));
    dblGreenRTSD_all = std(cell2mat(behave(lstGreenIndexes,4)));
    dblGreenRTSD_blk1 = std(cell2mat(behave(intersect(lstGreenIndexes,lstBlockIndexes1),4)));
    dblGreenRTSD_blk2 = std(cell2mat(behave(intersect(lstGreenIndexes,lstBlockIndexes2),4)));
    dblGreenRTSD_blk3 = std(cell2mat(behave(intersect(lstGreenIndexes,lstBlockIndexes3),4)));
    dblGreenRTSD_blk4 = std(cell2mat(behave(intersect(lstGreenIndexes,lstBlockIndexes4),4)));
    dblGreenRTSD_blk5 = std(cell2mat(behave(intersect(lstGreenIndexes,lstBlockIndexes5),4)));
    dblGreenRTSD_blk6 = std(cell2mat(behave(intersect(lstGreenIndexes,lstBlockIndexes6),4)));
    dblGreenRTSD_blk7 = std(cell2mat(behave(intersect(lstGreenIndexes,lstBlockIndexes7),4)));
    dblGreenRTSD_blk8 = std(cell2mat(behave(intersect(lstGreenIndexes,lstBlockIndexes8),4)));
    dblRedRTSD_all = std(cell2mat(behave(lstRedIndexes,4)));
    dblRedRTSD_blk1 = std(cell2mat(behave(intersect(lstRedIndexes,lstBlockIndexes1),4)));
    dblRedRTSD_blk2 = std(cell2mat(behave(intersect(lstRedIndexes,lstBlockIndexes2),4)));
    dblRedRTSD_blk3 = std(cell2mat(behave(intersect(lstRedIndexes,lstBlockIndexes3),4)));
    dblRedRTSD_blk4 = std(cell2mat(behave(intersect(lstRedIndexes,lstBlockIndexes4),4)));
    dblRedRTSD_blk5 = std(cell2mat(behave(intersect(lstRedIndexes,lstBlockIndexes5),4)));
    dblRedRTSD_blk6 = std(cell2mat(behave(intersect(lstRedIndexes,lstBlockIndexes6),4)));
    dblRedRTSD_blk7 = std(cell2mat(behave(intersect(lstRedIndexes,lstBlockIndexes7),4)));
    dblRedRTSD_blk8 = std(cell2mat(behave(intersect(lstRedIndexes,lstBlockIndexes8),4)));
    
    %find the subject entry (created during drug processing)
    idxSubjectEntry = find(str2double(strSubjectID) == matrixAllSubjectData(:,1));
    if (length(idxSubjectEntry) > 1)
        %more than 1 entry found, error
        error('More than 1 entry found in data matrix for subject %s',strSubjectID)
    elseif isempty(idxSubjectEntry)
        %subject not found
        fprintf('Warning: subject not found for %s, drug may not exist\n', strSubjectID)
        idxSubjectEntry = find(isnan(matrixAllSubjectData(:,1)),1);
        
        if isempty(idxNextEntry)
            %matrix is fully filled, so expand
            idxSubjectEntry = length(matrixAllSubjectData)+1;
        end
        
        matrixAllSubjectData(idxSubjectEntry,:) = NaN(1,size(matrixAllSubjectData,2));  %fill empty entries with NaNs instead of 0s
        matrixAllSubjectData(idxSubjectEntry,1) = str2double(strSubjectID);
    end
    
    %save into all subject data matrix
    matrixAllSubjectData(idxSubjectEntry,2) = dblGreenACC_blk1;
    matrixAllSubjectData(idxSubjectEntry,3) = dblGreenACC_blk2;
    matrixAllSubjectData(idxSubjectEntry,4) = dblGreenACC_blk3;
    matrixAllSubjectData(idxSubjectEntry,5) = dblGreenACC_blk4;
    matrixAllSubjectData(idxSubjectEntry,6) = dblGreenACC_blk5;
    matrixAllSubjectData(idxSubjectEntry,7) = dblGreenACC_blk6;
    matrixAllSubjectData(idxSubjectEntry,8) = dblGreenACC_blk7;
    matrixAllSubjectData(idxSubjectEntry,9) = dblGreenACC_blk8;
    matrixAllSubjectData(idxSubjectEntry,10) = dblGreenACC_all;
    matrixAllSubjectData(idxSubjectEntry,11) = dblRedACC_blk1;
    matrixAllSubjectData(idxSubjectEntry,12) = dblRedACC_blk2;
    matrixAllSubjectData(idxSubjectEntry,13) = dblRedACC_blk3;
    matrixAllSubjectData(idxSubjectEntry,14) = dblRedACC_blk4;
    matrixAllSubjectData(idxSubjectEntry,15) = dblRedACC_blk5;
    matrixAllSubjectData(idxSubjectEntry,16) = dblRedACC_blk6;
    matrixAllSubjectData(idxSubjectEntry,17) = dblRedACC_blk7;
    matrixAllSubjectData(idxSubjectEntry,18) = dblRedACC_blk8;
    matrixAllSubjectData(idxSubjectEntry,19) = dblRedACC_all;
    matrixAllSubjectData(idxSubjectEntry,38) = dblGreenRT_blk1;
    matrixAllSubjectData(idxSubjectEntry,39) = dblGreenRT_blk2;
    matrixAllSubjectData(idxSubjectEntry,40) = dblGreenRT_blk3;
    matrixAllSubjectData(idxSubjectEntry,41) = dblGreenRT_blk4;
    matrixAllSubjectData(idxSubjectEntry,42) = dblGreenRT_blk5;
    matrixAllSubjectData(idxSubjectEntry,43) = dblGreenRT_blk6;
    matrixAllSubjectData(idxSubjectEntry,44) = dblGreenRT_blk7;
    matrixAllSubjectData(idxSubjectEntry,45) = dblGreenRT_blk8;
    matrixAllSubjectData(idxSubjectEntry,46) = dblGreenRT_all;
    matrixAllSubjectData(idxSubjectEntry,47) = dblRedRT_blk1;
    matrixAllSubjectData(idxSubjectEntry,48) = dblRedRT_blk2;
    matrixAllSubjectData(idxSubjectEntry,49) = dblRedRT_blk3;
    matrixAllSubjectData(idxSubjectEntry,50) = dblRedRT_blk4;
    matrixAllSubjectData(idxSubjectEntry,51) = dblRedRT_blk5;
    matrixAllSubjectData(idxSubjectEntry,52) = dblRedRT_blk6;
    matrixAllSubjectData(idxSubjectEntry,53) = dblRedRT_blk7;
    matrixAllSubjectData(idxSubjectEntry,54) = dblRedRT_blk8;
    matrixAllSubjectData(idxSubjectEntry,55) = dblRedRT_all;
    matrixAllSubjectData(idxSubjectEntry,74) = dblGreenRTSD_blk1;
    matrixAllSubjectData(idxSubjectEntry,75) = dblGreenRTSD_blk2;
    matrixAllSubjectData(idxSubjectEntry,76) = dblGreenRTSD_blk3;
    matrixAllSubjectData(idxSubjectEntry,77) = dblGreenRTSD_blk4;
    matrixAllSubjectData(idxSubjectEntry,78) = dblGreenRTSD_blk5;
    matrixAllSubjectData(idxSubjectEntry,79) = dblGreenRTSD_blk6;
    matrixAllSubjectData(idxSubjectEntry,80) = dblGreenRTSD_blk7;
    matrixAllSubjectData(idxSubjectEntry,81) = dblGreenRTSD_blk8;
    matrixAllSubjectData(idxSubjectEntry,82) = dblGreenRTSD_all;
    matrixAllSubjectData(idxSubjectEntry,83) = dblRedRTSD_blk1;
    matrixAllSubjectData(idxSubjectEntry,84) = dblRedRTSD_blk2;
    matrixAllSubjectData(idxSubjectEntry,85) = dblRedRTSD_blk3;
    matrixAllSubjectData(idxSubjectEntry,86) = dblRedRTSD_blk4;
    matrixAllSubjectData(idxSubjectEntry,87) = dblRedRTSD_blk5;
    matrixAllSubjectData(idxSubjectEntry,88) = dblRedRTSD_blk6;
    matrixAllSubjectData(idxSubjectEntry,89) = dblRedRTSD_blk7;
    matrixAllSubjectData(idxSubjectEntry,90) = dblRedRTSD_blk8;
    matrixAllSubjectData(idxSubjectEntry,91) = dblRedRTSD_all;
    matrixAllSubjectData(idxSubjectEntry,110) = length(lstExcludedIndexes);
    matrixAllSubjectData(idxSubjectEntry,111) = length(lstExcludedIndexes)/size(behave,1)*100;
    
    clear behave structData
end
matrixAllSubjectData = sortrows(matrixAllSubjectData,1);    %infiles were not sorted correctly by ascending subject id

if flagDebug
    error('Debug complete')
end

%export to file

fprintf('\n\nSaving...')

fid = fopen([dirOutput,strOutfile],'w+');    %open file to write (overwrites)
fprintf(fid,'%s,',arrHeader{1:end-1});       %print headers
fprintf(fid,'%s\n',arrHeader{end});

%process all subject rows
for i = 1:size(matrixAllSubjectData,1)
    fprintf(fid,'%i,',matrixAllSubjectData(i,1));        %subject ID
    fprintf(fid,'%.3f,',matrixAllSubjectData(i,2:37));	 %accuracy
    fprintf(fid,'%.0f,',matrixAllSubjectData(i,38:73));	 %reaction time
    fprintf(fid,'%.3f,',matrixAllSubjectData(i,74:109)); %reaction time SD
    fprintf(fid,'%i,',matrixAllSubjectData(i,110));      %rejected trial count plac
    fprintf(fid,'%.1f,',matrixAllSubjectData(i,111));    %rejected trial percent plac
    fprintf(fid,'%i,',matrixAllSubjectData(i,112));      %rejected trial count drug
    fprintf(fid,'%.1f\n',matrixAllSubjectData(i,113));   %rejected trial percent drug
end

fclose(fid);
fprintf('done\n')

clear all