% Set Electrode Regions

%all PFC
all_PFC = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23];  
%PFC left
left_PFC = [1 4 6 7 8 15 16 17];
%PFC mid
mid_PFC = [2 9 10 11 18 19 20];
%PFC right
right_PFC = [3 5 12 13 14 21 22 23];
%all Central
all_Central = [24 25 26 27 28 29 30 31 32 34 35 36 37 38 39 40 41 42];
%Central left
left_Central = [24 25 26 34 35 36];
%Central mid
mid_Central = [27 28 29 37 38 39];
%Central right
right_Central = [30 31 32 40 41 42];
%all Posterior
all_Post = [44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 61 62 63];
%Posterior left
left_Post = [44 45 46 53 54];
%Posterior mid
mid_Post = [47 48 49 55 56 57 61 62 63];
%Posterior right
right_Post = [50 51 52 58 59];

%%
%dose_type = [{'controls'} {'patients'}];
dose_type = [{'controls'}];
%PFC_type = [{'all_PFC'} {'left_PFC'} {'mid_PFC'} {'right_PFC'}]
%PFC_type = [{'all_Central'} {'left_Central'} {'mid_Central'} {'right_Central'}]
%PFC_type = [{'all_Post'} {'left_Post'} {'mid_Post'} {'right_Post'}]
PFC_type = [{'left_PFC'} {'mid_PFC'} {'right_PFC'} {'left_Central'} {'mid_Central'} {'right_Central'} {'left_Post'} {'mid_Post'} {'right_Post'}]


time_interval = [1:1632]



for l = 1:length(PFC_type)
    
    eval(['PFC =', PFC_type{l}])

for k= 1:length(dose_type)

    dose = dose_type{k}
    
j=1;
dose_struct_all = eval(['ERSP.broadband.',dose]);
for i=1:length(dose_struct_all)
    if length(dose_struct_all{i}) ~= 0
        eval([dose,'_struct{j} = dose_struct_all{i}'])
        j = j+1
    end
end

clear i j dose_struct_all



%yield <subject X electrode X freqency X time_interval>

for i=1:length(eval([dose,'_struct']))
      for elec = 1:length(PFC)
        contrast_leads = eval([dose,'_struct{i}.order.lead']);
        %contrast_leads = eval([dose,'_struct{i}.nontarget.lead']);
        %contrast_leads = eval([dose,'_struct{i}.sub.lead']);
        if isempty(contrast_leads{PFC(elec)}) == 0
                theta_struct(i,elec,:,:) = (contrast_leads{PFC(elec)}(:,time_interval));
        end
        
    end
          
end

%squeeze the mean across electrode_type; yields <subject X frequency X
%time_interval>. squeeze again for mean <frequency X time_interval>

eval([PFC_type{l},'_',dose,'_theta_struct_initial = theta_struct;']);

mean_theta_struct = squeeze(nanmean(squeeze(nanmean(theta_struct,2))));
%mean_theta_struct = squeeze(nanmean(theta_struct,2));

eval([PFC_type{l},'_',dose,'_theta_struct = mean_theta_struct;']);


clearvars i k theta_struct mean_theta_struct
    
 
end

%eval([PFC_type{l},'_Contrast_theta_struct = ',PFC_type{l},'_controls_theta_struct - ',PFC_type{l},'_patients_theta_struct']);

clear PFC theta_struct Drug_struct Plac_struct
end

%% single subject

group =['CMET_broadband']

% % broadband
%  freqs = linspace(4,30,27)
%  time_interval =  linspace(-300,1300,161)
%'0gi
 
 % broadband 500 avg
%  freqs = linspace(4,30,27)
%  time_interval =  linspace(-220,2216,886)
 
 freqs = freqs_broadband
 time_interval =  times_broadband

theta_type = [{'controls'}]
%gamma_type = [{'Plac'}]
%PFC_type = [{'all_PFC'} {'left_PFC'} {'mid_PFC'} {'right_PFC'}]
%PFC_type = [{'all_Central'} {'left_Central'} {'mid_Central'} {'right_Central'}]
%PFC_type = [{'all_Post'} {'left_Post'} {'mid_Post'} {'right_Post'}]
PFC_type = [{'left_PFC'} {'mid_PFC'} {'right_PFC'} {'left_Central'} {'mid_Central'} {'right_Central'} {'left_Post'} {'mid_Post'} {'right_Post'}]


for k = 1:length(theta_type)
    
for j = 1:length(PFC_type)
    
    sjset1 = eval([PFC_type{j},'_',theta_type{k},'_theta_struct(1,:,:)']);
    sjset=squeeze(sjset1);

    
    
    close all

plot_name = [group,'_',theta_type{k},'_',PFC_type{j}]

%h = tftopo(sjset,time_interval,freqs,'title',[plot_name],'verbose','off');
h = tftopo(sjset,'logfreq','native');
    %caxis([-.33 .33])
    %caxis([-3 3])
    %caxis([-2 2])
    colorbar
    figure;
    %saveas(1, [plot_name,'.tif'],'tiff') %Use for publication
    saveas(1, [plot_name,'.jpg'],'jpg')
end
end
