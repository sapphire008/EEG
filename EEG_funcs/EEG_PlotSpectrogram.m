function EEG_PlotSpectrogram(data_mat,times_out,freqs_out,save_dir,...
    plot_name,varargin)
% Plot Spectrogram
% EEG_PlotSpectrogram(data_mat,times_out,freqs_out,save_dir,plot_name)
% Required Inputs: 
%       data_mat:   2-D data matrix
%       times_out:  time range
%       freqs_out:  frequency range
%       plot_name:  title of the plot
%
% Optional Inputs:
%       'ColorRange'(double):   colorbar range, imposing symmetry 
%                               (default:[-1.5,1.5])
%       'freqs_step' (double):  frequency increment (default:5)
%       'logfreqs' (char):      plot frequency using log scale
%                               (default: 'off')

%inspect variable inputs
flag = InspectVarargin(varargin,...
    {'ColorRange','freqs_step','logfreq','ReverseColor','ColorMap'},...
    {[-1.5,1.5],5,'off',0,[]});
%make sure save_dir ends with a forward slash
if ~strcmpi(save_dir(end),'/')
    save_dir = [save_dir,'/'];
end
%freqs x times plot
h = tftopo(data_mat,times_out,freqs_out,...
    'title',strrep(plot_name,'_','\_'),'verbose','off','cbar','on',...
    'logfreq',flag.logfreq);
caxis(flag.ColorRange);
set(gca,'YTick',[freqs_out(1):flag.freqs_step:freqs_out(end)]);
if ~isempty(flag.ColorMap)
    colormap(flag.ColorMap);
end
h1=colorbar;
if flag.ReverseColor
    colormap (flag.ColorMap(end:-1:1,:));
    set(h1,'YDir','reverse');
end

saveas(gcf, [save_dir,plot_name,'.tif'],'tiff');
close all;

end