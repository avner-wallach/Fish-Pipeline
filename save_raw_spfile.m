function []=save_raw_spfile(amp,changroups,eodind,b_dhpf,a_dhpf,sesspath)
%   GET_SPIKE_TIMES extract spike times from continuous tetrode recording.
%   t- input time vector
%   amp- input voltage matrix (column/channel)
%   changroups- cell array with vectors of channel ids for each tetrode to
%               be analyzed.
%   eodtimes- times of eod pulses
%% PARAMETERS
samplerate=str2num(getenv('SAMPLERATE'));
eodtimes=eodind/samplerate;

%blanking
blank_gap=str2num(getenv('BLANKGAP'))*samplerate;
blank_pre=str2num(getenv('BLANKPRE'))*samplerate;
blank_post=str2num(getenv('BLANKPOST'))*samplerate;

%spike detection
reftime=str2num(getenv('REFTIME'));
spikewidth=str2num(getenv('SPIKEWIDTH'));
minpoint=str2num(getenv('MINPOINT'));
thfactor=str2num(getenv('THFACTOR'));
artth=str2num(getenv('ARTTH'));
adcthreshold=str2num(getenv('ADCTHRESHOLD'));
direction=getenv('DIRECTION');

%clustering
blocksize0=str2num(getenv('SBLOCKSIZE')); %min clustering block size
clustnum=str2num(getenv('CLUSTNUM')); %number of clusters in each block
kmeansdims=str2num(getenv('KMEANSDIMS')); %number of clusters in each block
isiedges=[0:.1:5];
isibins=edge2bin(isiedges);

%% 
K=10;
L=10;
I=eodind*ones(1,blank_pre+blank_post+1) + ones(size(eodind,1),1)*[-blank_pre:blank_post];
Io=eodind*ones(1,K*2) + ones(size(eodind,1),1)*[[(-blank_pre-K*L):L:(-blank_pre-1)] [(blank_post+1):L:(blank_post+K*L)]];
I_gap=eodind*ones(1,blank_pre+blank_post+1+2*blank_gap) + ones(size(eodind,1),1)*[-(blank_pre+blank_gap):(blank_post+blank_gap)];
I=unique(I(:));
Io=unique(Io(:));
I_gap=unique(I_gap(I_gap>0 & I_gap<=size(amp,1)));
Io=setdiff(Io,intersect(I,Io));
I=I(I>0 & I<size(amp,1));
Io=Io(Io>0 & Io<size(amp,1));

C=zeros(size(amp,1),numel(cell2mat(changroups)));
k=0;
for i=1:numel(changroups)
    A=amp(:,changroups{i}); %take channels of current tetrode    
    for j=1:numel(changroups{i})
        a=A(:,j);
        %blank EODs
        B=interp1(Io(:),a(Io(:)),I,'pchip');
        a(I)=B;
        a=FiltFiltM(b_dhpf,a_dhpf,a);
        %blank EODs +gaps
        a(I_gap)=0;
        A(:,j)=a;
    end    
    C(:,k+[1:size(A,2)])=A;
    k=k+size(A,2);
end

%append to binary file
F=fopen(filename,'a');
fwrite(F,C','int16');
fclose(F);
end