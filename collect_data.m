function[frame,eod,file]=collect_data(segnums,frame,eod,file)
%COLLECT_DATA collect data from segments numbered segnums, taking
%physiology data from channels chnind

datapath=getenv('DATAPATH');
sdate=getenv('SESSDATE');
samplerate=str2num(getenv('SAMPLERATE'));
trace_b=str2num(getenv('TRACEBLNK')); %trace blank time, ms
% outchans=str2num(getenv('OUTCHANS')); %vector of channels to keep

trwind=str2num(getenv('TRACEWIND')); %window for measuring field
tracedown=str2num(getenv('TRACEDOWN')); %direction of LFP 1=down, 0=up
tracedlg=getenv('TRACEDLG'); %trace dialog
% trwind=[1 5];
% tracedown=1;
traceth=str2num(getenv('TRACETH')); %threshold for cleaning noise
tracenum=str2num(getenv('TRACENUM'));    %number of averaged traces
tracevar=getenv('TRACEVAR');   %ordering variable for averaged traces

dtrend=getenv('DETREND');
dtrendt=str2num(getenv('DETRENDT'));
glitchth=str2num(getenv('GLITCHTH'));
glitcht=str2num(getenv('GLITCHT'));
glitchref=str2num(getenv('GLITCHREF'));

sesspath=[datapath,'\',sdate,'\'];
if(nargin==1)
    frame.t=[];
    frame.data=[];
    frame.ind0=[1];
    eod.t=[];
    eod.data=[];
    eod.ind0=[1];
end
    
latnames={};
ampnames={};
for i=1:numel(segnums)
    i
    if(~exist([sesspath,'data_',num2str(segnums(i)),'.mat']))
        continue;
    end
    s=load([sesspath,'data_',num2str(segnums(i))]);
    data=s.data;
    N=min([numel(data.FRAME.t),size(data.FRAME.posture,1),numel(data.FRAME.accx),numel(data.FRAME.accy),numel(data.FRAME.accz)]);
    frame.t=[frame.t;data.FRAME.t(1:N)];
    frame.data=[frame.data;data.FRAME.posture(1:N,:) data.FRAME.accx(1:N) data.FRAME.accy(1:N) data.FRAME.accz(1:N)];
    frame.fnames={data.FILE.model{:} 'accx' 'accy' 'accz'};
    frame.ind0=[frame.ind0;size(frame.data,1)];
    if(i==1)
%         F=plot_traces;   
%         if(~numel(outchans))
%             [chnind]=chndlg()
%         else
        outchans=[1:size(data.EOD.traces,3)];
           chnind=outchans;
%         end
%         close(F);
        for k=1:numel(chnind)
            latnames{k}=['lat',num2str(chnind(k))];
            ampnames{k}=['amp',num2str(chnind(k))];
            aucnames{k}=['auc',num2str(chnind(k))];
%             crnames{k}=['cr',num2str(chnind(k))];        
        end        
        if(numel(trwind)>2)
%             for j=1:numel(chnind)
%                 F=plot_traces(chnind(j));
%                 [win(j,:),pol(j)]=lfpdlg(j);
%                 close(F);                           
%             end            
            for j=1:numel(chnind)
                win(j,:)=trwind(j,:);
                pol(j)=(~tracedown(j));
            end

        else        
            for j=1:numel(chnind)
                win(j,:)=trwind;
                pol(j)=(~tracedown);
            end
        end
%             lfp_temp{j}=get_lfp_shape(data.EOD.traces(:,:,chnind(j)),win(j,:));

        if(nargin==1)
            file=data.FILE;
            file.offset=[];
            file.title=[];
        end
    end
    [lat,amp,auc,avtraces(:,:,:,i)]=get_lfp_stats(data.EOD.traces);
    % get median readings (standatized) from entire site
    
    gind=get_good_inds(data.EOD.t,data.EOD.traces,amp);
    lat=lat(gind,:);
    amp=amp(gind,:);
    auc=auc(gind,:);
%     all_lat=nanmedian(zscore(lat(:,[2 4]),[],1),2);
%     all_amp=nanmedian(zscore(amp(:,[2 4]),[],1),2);

    deodt=[nan diff(data.EOD.t)];
    eod.t=[eod.t;data.EOD.t(gind)'];
    iei=deodt(gind)';    
    eod.data=[eod.data;iei data.EOD.posture(gind,:) data.EOD.accx(gind)' data.EOD.accy(gind)' data.EOD.accz(gind)' lat amp auc];% cr(gind,:)];
    eod.fnames={'iei' data.FILE.model{:} 'accx' 'accy' 'accz' latnames{:} ampnames{:} aucnames{:}};%crnames{:}};
    eod.ind0=[eod.ind0;size(eod.data,1)];
    
    file.offset=[file.offset;data.FILE.offset];    
    
end

eod.avtraces=avtraces;

function F=plot_traces(ind)
    K=floor(sqrt(numel(outchans)));
    M=ceil(numel(outchans)/K);
    F=figure;
    x=[1:size(data.EOD.traces,2)]/samplerate*1e3+trace_b;
    if(nargin>0)
        q=quantile(data.EOD.traces(:,:,ind),[0.25 0.5 0.75]);
        my_plotWithConf_Q(x,q',[0 0 0],1);
    else
        for j=1:numel(outchans)
            q=quantile(data.EOD.traces(:,:,outchans(j)),[0.25 0.5 0.75]);
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
    if(tracedown)
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

function [lat,amp,auc,avtraces]=get_lfp_stats(traces)
    x=[1:size(traces,2)]/samplerate*1e3+trace_b;    
    t_align=dtrendt-trace_b;
    N_align=ceil(t_align*samplerate/1e3);
%     N_align=[12];

    %lpf filter
    lp=100*2/samplerate;
    [NN, Wn] = buttord( lp, lp .* [0.75], 3, 20);
    [B1,A1] = butter(NN,Wn,'high');
    
    for k=1:numel(chnind)
        
        if(numel(glitchth))        
            t_noise=glitcht;
            N_noise=min(max(ceil(t_noise*samplerate/1e3),1),size(traces,2));            
            ind=find(max(abs(diff(traces(:,N_noise(1):N_noise(2),chnind(k)),1,2)),[],2)>glitchth(k));
            if(numel(ind))
                traces(ind,:,chnind(k))=nan;
            end
        end

        %clean glitches
        %detrend        
        if(strcmp(dtrend,'on'))        
            traces(:,:,chnind(k))=my_detrend(traces(:,:,chnind(k)),N_align);
        elseif(strcmp(dtrend,'matlab'))        
            if(numel(N_align))
                traces(:,:,chnind(k))=detrend(traces(:,:,chnind(k))','linear',N_align)';
            else
                traces(:,:,chnind(k))=detrend(traces(:,:,chnind(k))')';
            end
        elseif(strcmp(dtrend,'hpf'))        
            traces(:,:,chnind(k))=filtfilt(B1,A1,traces(:,:,chnind(k))')';
        end
        
        ind=find(x>=win(k,1) & x<=win(k,2));
        [M,Mind]=max(traces(:,ind,chnind(k)),[],2);
        [m,mind]=min(traces(:,ind,chnind(k)),[],2);
        if(pol(k))
            amp(:,k)=M;
            auc(:,k)=sum(traces(:,ind,chnind(k)).*(traces(:,ind,chnind(k))>0),2);
            lind=Mind;
        else
            amp(:,k)=-m;
            auc(:,k)=-sum(traces(:,ind,chnind(k)).*(traces(:,ind,chnind(k))<0),2);
            lind=mind;            
        end        
%         amp(:,k)=M-m;
%         lind=min(Mind,mind);
        lat(:,k)=lind/samplerate*1e3+win(k,1);
        I=find(lind==1 | lind==numel(ind)); %no local extremum
        amp(I,k)=NaN;
        lat(I,k)=NaN;
        
        if(strcmp(tracevar,'amp'))
            [x,avtraces(:,:,k)]=plot_graded(x,traces(:,:,chnind(k)),amp(:,k),tracenum,0);
        elseif(strcmp(tracevar,'lat'))
            [x,avtraces(:,:,k)]=plot_graded(x,traces(:,:,chnind(k)),lat(:,k),tracenum,0);
        else
            [x,avtraces(:,:,k)]=plot_graded(x,traces(:,:,chnind(k)),auc(:,k),tracenum,0);
        end

    end
end

function indg=get_good_inds(t,traces,amps)
%     indg=find(~isnan(sum(sum(traces,3),2)));    %good traces
%     indb=find(isnan(sum(sum(traces,3),2)));     %bad traces
    A=prod(prod(isnan(traces),3),2);
    indg=find(~A);    %good traces
    indb=find(A);     %bad traces
    %remove within refractory period after noise
    for j=1:numel(indb)
        T=(t(indg)-t(indb(j)));
        idx=find(T>0 & T<glitchref);
        indg(idx)=[];
    end
    
    %remove no signal traces
    indg=intersect(indg,find(prod(abs(amps)>traceth,2)));
end
    
end