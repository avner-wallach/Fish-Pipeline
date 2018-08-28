function [traces,tr_amp,tr_lat]=clean_traces(traces)
%% get params
trace_b=str2num(getenv('TRACEBLNK')); %trace blank time, ms
trace_T=str2num(getenv('TRACET')); %total trace duration
samplerate=str2num(getenv('SAMPLERATE'));
T=[0:size(traces,2)-1]/samplerate*1e3;
sthreshold=str2num(getenv('TRACETH')); %std threshold for cleaning noise
% sthreshold=1000;
trwind=str2num(getenv('TRACEWIND')); %window for measuring field
troffwind=str2num(getenv('TRACEOFFWIND')); %window for removing offsets
tracedown=str2num(getenv('TRACEDOWN')); %direction of LFP 1=down, 0=up

%% remove offset
if(numel(troffwind))
    T0=troffwind(1);
    T1=troffwind(2);
    ind=find(T>T0 & T<T1);
    t0=mean(traces(:,ind,:),2);
    traces=traces-repmat(t0,1,size(traces,2),1);
end
%% trace window
T0=trwind(1);
T1=trwind(2);
ind=find(T>T0 & T<T1);

%% remove noise
if(numel(sthreshold))
    m=std(traces(:,:,:),[],2);
    for i=1:size(traces,3)
        idx=find(m(:,:,i)>sthreshold);
        traces(idx,:,i)=nan;
    end
end
%% find min
[m,i]=min(traces(:,ind,:),[],2);
[M,I]=max(traces(:,ind,:),[],2);
M(:,:,tracedown==1)=m(:,:,tracedown==1);
I(:,:,tracedown==1)=i(:,:,tracedown==1);
tt=T(ind);
tr_amp=permute(M,[1 3 2]);
tr_lat=permute(tt(I),[1 3 2]);
end