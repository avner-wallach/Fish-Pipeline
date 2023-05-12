function []=plot_raw_traces(amp,ops,segnum,eodind,t)
%   GET_SPIKE_TIMES extract spike times from continuous tetrode recording.
%   t- input time vector
%   amp- input voltage matrix (column/channel)
%   changroups- cell array with vectors of channel ids for each tetrode to
%               be analyzed.
%   eodtimes- times of eod pulses
%% PARAMETERS
samplerate=ops.samplerate;

%EOD blanking
eod_pre=ops.eodblankpre*samplerate;
eod_post=ops.eodblankpost*samplerate;

%LFP blanking
blank_pre=ops.blankpre*samplerate;
blank_post=ops.blankpost*samplerate;

%changroups
changroups=cell2mat(ops.seg(segnum).LFPgroups);
A=amp(:,changroups); %take channels of current tetrode    

%% new method filtering : median filtering + EOD blanking + LFP template cancelation
K=2e-3*samplerate; %med filt kernel
A=A-medfilt1(A,K,[],1); %use median filter to remove slow components

%EOD blanking sections
I_eod=eodind*ones(1,eod_pre+eod_post+1) + ones(size(eodind,1),1)*[-eod_pre:eod_post];
I_eod=unique(I_eod(:));
I_eod=I_eod(I_eod>0 & I_eod<size(amp,1));

A(I_eod(:),:)=0; %blank eod

%LFP blanking sections
eodind=eodind(eodind>blank_pre & eodind<(size(amp,1)-blank_post));
I=eodind*ones(1,blank_pre+blank_post+1) + ones(size(eodind,1),1)*[-blank_pre:blank_post];

for i=1:size(A,2)
    a=A(:,i);
    Y=a(I)';    
    [x,tr]=plot_graded(1:size(I,2),Y',[],max(abs(Y)),50,0);
    cc=corr(Y,tr');
    [c0,mm]=max(cc,[],2);
    T=tr(mm,:)';
    beta=(sum(T.*Y))./vecnorm(T).^2;
    D=beta'*ones(1,size(tr,2)).*T';    
    a(I)=a(I)-D;
    A(:,i)=a;     
end

%% plot data
A(1,:)=0;
% A=A/max(abs(A(:)));
A=A/100;
for i=1:numel(ops.seg(segnum).LFPgroups)
    plot(t,A(:,ops.seg(segnum).LFPgroups{i})+i);
    s=std(A(:,ops.seg(segnum).LFPgroups{i}),[],1);    
    hold on;
    H=plot([t(1) t(end)],(i+max(s)*3)*[1 1],[t(1) t(end)],(i-max(s)*3)*[1 1]);
    set(H,'Color',0.5*[1 1 1]);
end
end