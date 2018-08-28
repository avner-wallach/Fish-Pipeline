function[frame,eod,file]=collect_data(segnums)
%COLLECT_DATA collect data from segments numbered segnums, taking
%physiology data from channels chnind

datapath=getenv('DATAPATH');
sdate=getenv('SESSDATE');
samplerate=str2num(getenv('SAMPLERATE'));
trace_b=str2num(getenv('TRACEBLNK')); %trace blank time, ms
outchans=str2num(getenv('OUTCHANS')); %vector of channels to keep
dtrendt=str2num(getenv('DETRENDT'));

trwind=str2num(getenv('TRACEWIND')); %window for measuring field
tracedown=str2num(getenv('TRACEDOWN')); %direction of LFP 1=down, 0=up
traceth=str2num(getenv('TRACETH')); %threshold for cleaning noise

sesspath=[datapath,'\',sdate,'\'];
frame.t=[];
frame.data=[];
frame.ind0=[1];
eod.t=[];
eod.data=[];
eod.ind0=[1];
latnames={};
ampnames={};
for i=1:numel(segnums)
    i
    s=load([sesspath,'data_',num2str(segnums(i))]);
    data=s.data;
    N=min([numel(data.FRAME.t),size(data.FRAME.posture,1),numel(data.FRAME.roll),numel(data.FRAME.pitch),numel(data.FRAME.yaw)]);
    frame.t=[frame.t;data.FRAME.t(1:N)];
    frame.data=[frame.data;data.FRAME.posture(1:N,:) data.FRAME.roll(1:N) data.FRAME.pitch(1:N) data.FRAME.yaw(1:N)];
    frame.fnames={data.FILE.model{:} 'roll' 'pitch' 'yaw'};
    frame.ind0=[frame.ind0;size(frame.data,1)];
    if(i==1)
        F=plot_traces;   
%         [chnind]=chndlg()
        chnind=outchans;
        close(F);
        for k=1:numel(chnind)
            latnames{k}=['lat',num2str(chnind(k))];
            ampnames{k}=['amp',num2str(chnind(k))];
%             crnames{k}=['cr',num2str(chnind(k))];        
        end        
        for j=1:numel(chnind)
            F=plot_traces(chnind(j));
%             [win(j,:),pol(j)]=lfpdlg(j);
            win(j,:)=trwind(j,:);
            pol(j)=(tracedown(j)=='+');

            lfp_temp{j}=get_lfp_shape(data.EOD.traces(:,:,chnind(j)),win(j,:));
            close(F);
        end
        file=data.FILE;
    end
    [lat,amp,cr]=get_lfp_stats(data.EOD.traces);
    % get median readings (standatized) from entire site
    
    gind=get_good_inds(data.EOD.t,data.EOD.traces,amp);
    lat=lat(gind,:);
    amp=amp(gind,:);
    all_lat=nanmedian(zscore(lat(:,[2 4]),[],1),2);
    all_amp=nanmedian(zscore(amp(:,[2 4]),[],1),2);

    deodt=[nan diff(data.EOD.t)];
    eod.t=[eod.t;data.EOD.t(gind)'];
    iei=deodt(gind)';    
    eod.data=[eod.data;iei data.EOD.posture(gind,:) data.EOD.roll(gind)' data.EOD.pitch(gind)' data.EOD.yaw(gind)' lat amp all_lat all_amp];% cr(gind,:)];
    eod.fnames={'iei' data.FILE.model{:} 'roll' 'pitch' 'yaw' latnames{:} ampnames{:} 'all_lat' 'all_amp'};%crnames{:}};
    eod.ind0=[eod.ind0;size(eod.data,1)];
    
    file.offset(i)=data.FILE.offset;    
    
end

function F=plot_traces(ind)
    K=floor(sqrt(numel(outchans)));
    M=numel(outchans)/K;
    F=figure;
    x=[1:size(data.EOD.traces,2)]/samplerate*1e3+trace_b;
    if(nargin>0)
        q=quantile(data.EOD.traces(:,:,ind),[0.25 0.5 0.75]);
        my_plotWithConf_Q(x,q',[0 0 0],1);
    else
        for j=1:numel(outchans)
            q=quantile(data.EOD.traces(:,:,j),[0.25 0.5 0.75]);
            subplot(K,M,j), my_plotWithConf_Q(x,q',[0 0 0],1);
        end
    end
end

function [chnind]=chndlg()
    prompt = {'chan numbers:'};
    dlg_title = '';
    num_lines = 1;
    defaultans = {'1 2 3 4'};
    answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
    chnind=str2num(answer{1});
end

function [win,pol]=lfpdlg(idx)
    prompt = {'window start:','window end:','polarity:'};
    dlg_title = 'LFP window params';
    num_lines = 1;
    if(tracedown(idx))
        trdir='-';
    else
        trdir='+';
    end
    defaultans = {num2str(trwind(1)) num2str(trwind(2)) trdir};
    answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
    win=[str2num(answer{1}) str2num(answer{2})];
    pol=(answer{3}=='+');
end

function lfp_tmp=get_lfp_shape(traces,wind)
    x=[1:size(traces,2)]/samplerate*1e3+trace_b;        
    ind=find(x>=wind(1) & x<=wind(2));
    lfp_tmp=nanmedian(traces(:,ind),1);    
end

function [lat,amp,cr]=get_lfp_stats(traces)
    x=[1:size(traces,2)]/samplerate*1e3+trace_b;    
    t_align=dtrendt-trace_b;
    N_align=ceil(t_align*samplerate/1e3);
%     N_align=[12];
    
    for k=1:numel(chnind)
        %detrend
        traces(:,:,chnind(k))=detrend(traces(:,:,chnind(k))','linear',N_align)';
        
        ind=find(x>=win(k,1) & x<=win(k,2));
        [M,Mind]=max(traces(:,ind,chnind(k)),[],2);
        [m,mind]=min(traces(:,ind,chnind(k)),[],2);
        if(pol(k))
            amp(:,k)=M;
            lind=Mind;
        else
            amp(:,k)=-m;
            lind=mind;            
        end        
%         amp(:,k)=M-m;
%         lind=min(Mind,mind);
        lat(:,k)=lind/samplerate*1e3+win(k,1);
        I=find(lind==1 | lind==numel(ind)); %no local extremum
        amp(I,k)=NaN;
        lat(I,k)=NaN;
        
        %correlation with template
%         for m=1:size(traces,1)
%             cr(m,k)=corr(traces(m,ind,chnind(k))',lfp_temp{k}');
%         end
        cr=[];

    end
end

function indg=get_good_inds(t,traces,amps)
    indg=find(~isnan(sum(sum(traces,3),2)));    %good traces
    indb=find(isnan(sum(sum(traces,3),2)));     %bad traces
    %remove within refractory period after noise
    for j=1:numel(indb)
        T=(t(indg)-t(indb(j)));
        idx=find(T>0 & T<.15);
        indg(idx)=[];
    end
    
    %remove no signal traces
    indg=intersect(indg,find(prod(abs(amps)>traceth,2)));
end
    
end