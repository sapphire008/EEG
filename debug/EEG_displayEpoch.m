function EEG_displayEpoch(EEG)
for ep = 1:length(EEG.epoch)
    disp(['Epoch #',ep,EEG.epoch(ep).eventtype]);
end
end