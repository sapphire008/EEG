function avg_img=EEG_meanSpect(Spectogram_cell)
%Spectogram_cell must be a cell array of spectograms in a row
num_img=length(Spectogram_cell);%nubmer of images
%place holding an average image
avg_img=zeros(size(Spectogram_cell{1,1}));

% Calculate average image
for img=1:num_img
    avg_img=avg_img+Spectogram_cell{1,img};
end
%return the average image
disp(['Averaging across ', num2str(num_img), ' images']);
avg_img=avg_img./num_img;
end