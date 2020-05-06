function [spikestruct]=get_spike_times(amp,changroups,eodind,b_dhpf,a_dhpf,offset,adc,sesspath,daynum)
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
    
    if(strcmp(direction,'falling'))
        [spiketimes,spiketraces,spikeidx]=MySpikeFinder(min(A,[],2),samplerate,reftime,spikewidth,minpoint,thfactor,'falling',0);
    elseif(strcmp(direction,'rising'))
        [spiketimes,spiketraces,spikeidx]=MySpikeFinder(max(A,[],2),samplerate,reftime,spikewidth,minpoint,thfactor,'rising',0);
    else %both
        [spiketimes,spiketraces,spikeidx]=MySpikeFinder(max(abs(A),[],2),samplerate,reftime,spikewidth,minpoint,thfactor,'falling',0);
    end
    
    %adc is used as noise control
    if(exist('adc') & numel(adc))
        adctraces=adc(spikeidx)-mean(adc);
        indclean=find(max(abs(adctraces),[],1)<adcthreshold);
        spiketimes=spiketimes(indclean);
        spikeidx=spikeidx(:,indclean);
    end
    
    %get blocksize
    k=1;
    blocksize=numel(spiketimes);
    while((blocksize)>blocksize0)
        k=k+1;
        blocksize=numel(spiketimes)/k;        
    end
    
    %partition into blocks, cluster and save cluster stats     
    i1=1; i2=floor(blocksize);
    k=1;
    stimes=[];
    while(i1<numel(spiketimes))
        idx=i1:i2;
        sptraces=zeros(numel(idx),size(spikeidx,1)*size(A,2));
        for n=1:size(A,2)
            I=[1:size(spikeidx,1)]+(n-1)*size(spikeidx,1);
            a=A(:,n);
            sptraces(:,I)=a(spikeidx(:,idx)');
        end

        %remove artifacts   
        M=max(abs(sptraces),[],2);
        notart=find(M<artth);
        sptraces=sptraces(notart,:);
        idx=idx(notart);
        
        stimes=[stimes;spiketimes(idx)'];

        %pca
        [score]=pca(sptraces');
        
        %cluster to clustnum clusters
        [spikestruct(i).blocks(k).indices]=kmeans(score(:,kmeansdims),clustnum);
        
        for j=1:clustnum
            cidx=find(spikestruct(i).blocks(k).indices==j);
            %centroid            
            spikestruct(i).blocks(k).centroids(j,:)=mean(sptraces(cidx,:),1);
            isi=diff(spiketimes(idx(cidx)))*1e3;
            spikestruct(i).blocks(k).isihist(j,:)=histcounts(isi,isiedges,'Normalization','probability');
        end
        i1=i2+1;    k=k+1;
        i2=floor(k*blocksize);
    end
            
    spikestruct(i).raster=stimes-offset;
    spikestruct(i).isibins=isibins;
    
end

if(numel(sesspath))
    filename=[sesspath,'spdata',num2str(daynum)];
    %append to binary file
    F=fopen(filename,'a');
    fwrite(F,C','int16');
    fclose(F);
end

end