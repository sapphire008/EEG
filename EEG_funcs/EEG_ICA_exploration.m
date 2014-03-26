function EEG_ICA_exploration(EEG)
%input EEG data after ICA
    headplot('setup',EEG.chanlocs,'headplot.spl');
    headplot(EEG.data,'headplot.spl','view',[0 90],...
        'maplimits','absmax','cbar',0,'electrodes','off');
end