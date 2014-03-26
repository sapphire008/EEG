
datadir='/nfs/to-eeg/';  
subno=1;
addpath('/nfs/to-eeg/code/functions/')

for n=1:length(EEG.epoch)
Event_Order{n} = EEG.epoch(n).eventtype{1}; 
end

item_hit=strsearch('200',Event_Order);
order_hit=strsearch('220',Event_Order);

EEG=pop_select(EEG, 'nochannel', {'CB1' 'CB2' 'HEO' 'VEO'}); 
electrode_locations = cell2mat({EEG.chanlocs.type});

for n=1:length(electrode_locations)
    elec = electrode_locations(n);

 disp(['subject = ', num2str(subno)]) 
 disp(['electrode = ', num2str(elec)]) 
 times_array = [-150:10:17000];
  
          [tmp_ersp itc powbase times_broadband freqs_broadband]=newtimef( {EEG.data(n,:,order_hit) EEG.data(n,:,item_hit)}, ...
              EEG.pnts,[EEG.xmin EEG.xmax]*1000, EEG.srate,[6] ,'timesout', times_array,'freqs',[4 30], ...
              'nfreqs',27, 'elocs', EEG.chanlocs, 'baseline', [-150 0],'chaninfo', EEG.chaninfo,...
              'plotphase','off','plotersp','off','plotitc','off'); 
          
%           [tmp_ersp itc powbase times_broadband freqs_broadband]=newtimef( {EEG.data(n,:,order_hit) EEG.data(n,:,item_hit)}, ...
%               EEG.pnts,[EEG.xmin EEG.xmax]*1000, EEG.srate,[6] ,'timesout', times_array,'freqs',[4 30],'freqscale','log', ...
%               'nfreqs',27, 'elocs', EEG.chanlocs, 'baseline', [-150 0],'chaninfo', EEG.chaninfo,...
%               'plotphase','off','plotersp','off','plotitc','off'); 
                    
         order_ersp = tmp_ersp{1,1};
         item_ersp   = tmp_ersp{1,2};
         order_vs_item_ersp    = tmp_ersp{1,3};

         order_itc = itc{1,1};
         item_itc   = itc{1,2};
         order_vs_item_itc    = itc{1,3};

    ERSP.broadband.controls{subno}.order.lead{elec}=order_ersp;
    ERSP.broadband.controls{subno}.item.lead{elec}=item_ersp;
    ERSP.broadband.controls{subno}.order_vs_item.lead{elec}=order_vs_item_ersp;  
    ERSP.ITC.broadband.controls{subno}.order.lead{elec}=order_itc;
    ERSP.ITC.broadband.controls{subno}.item.lead{elec}=item_itc;
    ERSP.ITC.broadband.controls{subno}.order_vs_item.lead{elec}=order_vs_item_itc;
    
end

 

cd (datadir)
cd analysis

save('ERSP_broadband_toc001_no_log.mat','ERSP','times_broadband','freqs_broadband','-v7.3')





