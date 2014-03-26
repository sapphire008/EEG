c=1;
for n = 1:length(EEG.urevent)
  events(n)=EEG.urevent(n).type;
  latency(n)=EEG.urevent(n).latency;
end
events_mat=reshape(events',3,length(EEG.urevent)/3)';
latency_mat=reshape(latency',3,length(EEG.urevent)/3)';
latency_mat=diff(latency_mat,1,2);

%%
for n = 1:length(EEG.epoch)
    A(n)=length(EEG.epoch(n).event);
    if A~=2;
        disp(A(n));
    end
end