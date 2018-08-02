function[frame,eod,file]=collect_data( segnums)
%COLLECT_DATA collect data from segments numbered segnums, taking
%physiology data from channels chnind

datapath=getenv('DATAPATH');
sdate=getenv('SESSDATE');
samplerate=str2num(getenv('SAMPLERATE'));
trace_b=str2num(getenv('TRACEBLNK')); %trace blank time, ms

sesspath=[datapath,'\',sdate,'\'];
frame.t=[];
frame.data=[];
frame.ind0=[1];
eod.t=[];
eod.data=[];
eod.ind0=[1];
% chnind=7;
latnames={};
ampnames={};
for i=1:numel(segnums)
    i
    s=load([sesspath,'data_',num2str(segnums(i))]);
    data=s.data;
    
    frame.t=[frame.t;data.FRAME.t];
    frame.data=[frame.data;data.FRAME.posture data.FRAME.roll data.FRAME.pitch data.FRAME.yaw];
    frame.fnames={data.FILE.model{:} 'roll' 'pitch' 'yaw'};
    frame.ind0=[frame.ind0;size(frame.data,1)];
    if(i==1)
        F=plot_traces;   
        [chnind]=chndlg()
        close(F);
        for k=1:numel(chnind)
            latnames{k}=['lat',num2str(chnind(k))];
            ampnames{k}=['amp',num2str(chnind(k))];
            crnames{k}=['cr',num2str(chnind(k))];
        end        
        for j=1:numel(chnind)
            F=plot_traces(chnind(j));
            [win(j,:),pol(j)]=lfpdlg();
            lfp_temp{j}=get_lfp_shape(data.EOD.traces(:,:,chnind(j)),win(j,:))
            close(F);
        end
        file=data.FILE;
    end
    [lat,amp,cr]=get_lfp_stats(data.EOD.traces);

    eod.t=[eod.t;data.EOD.t'];
    iei=[NaN diff(data.EOD.t)];
    eod.data=[eod.data;iei' data.EOD.posture data.EOD.roll' data.EOD.pitch' data.EOD.yaw' lat amp cr];
    eod.fnames={'iei' data.FILE.model{:} 'roll' 'pitch' 'yaw' latnames{:} ampnames{:} crnames{:}};
    eod.ind0=[eod.ind0;size(eod.data,1)];
    
    file.offset(i)=data.FILE.offset;    
    
end

%pca
X=frame.data(:,4:13);
Xm=nanmean(X,1);
[coefs,score,latent,tsq,explained]=pca(X);
for i=1:size(coefs,1)
    pcanames{i}=['pca',num2str(i)];
end
frame.data=[frame.data score];
frame.fnames={frame.fnames{:} pcanames{:}};
file.pcacoefs=coefs;
file.explained=explained;
%EOD data pca
E=eod.data(:,5:14);
E=E-ones(size(E,1),1)*Xm; %center
Escore=E*coefs;
eod.data=[eod.data Escore];
eod.fnames={eod.fnames{:} pcanames{:}};

function F=plot_traces(ind)
    F=figure;
    x=[1:size(data.EOD.traces,2)]/samplerate*1e3+trace_b;
    if(nargin>0)
        q=quantile(data.EOD.traces(:,:,ind),[0.25 0.5 0.75]);
        my_plotWithConf_Q(x,q',[0 0 0],1);
    else
        for j=1:16
            q=quantile(data.EOD.traces(:,:,j),[0.25 0.5 0.75]);
            subplot(4,4,j), my_plotWithConf_Q(x,q',[0 0 0],1);
        end
    end
end

function [chnind]=chndlg()
    prompt = {'chan numbers:'};
    dlg_title = '';
    num_lines = 1;
    defaultans = {'5 7 9 14'};
    answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
    chnind=str2num(answer{1});
end

function [win,pol]=lfpdlg()
    prompt = {'window start:','window end:','polarity:'};
    dlg_title = 'LFP window params';
    num_lines = 1;
    defaultans = {'2','4','-'};
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
    for k=1:numel(chnind)
        ind=find(x>=win(k,1) & x<=win(k,2));
        if(pol(k))
            [amp(:,k),lind]=max(traces(:,ind,chnind(k)),[],2);
        else
            [amp(:,k),lind]=min(traces(:,ind,chnind(k)),[],2);
        end        
        lat(:,k)=lind/samplerate*1e3+trace_b;
        I=find(lind==1 | lind==numel(ind)); %no local extremum
        amp(I,k)=NaN;
        lat(I,k)=NaN;
        
        %correlation with template
        for m=1:size(traces,1)
            cr(m,k)=corr(traces(m,ind,chnind(k))',lfp_temp{k}');
        end

    end
end
end