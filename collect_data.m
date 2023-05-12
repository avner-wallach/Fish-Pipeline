function[frame,eod,file]=collect_data(segnums,frame,eod,file)
%COLLECT_DATA collect data from segments numbered segnums, taking
%physiology data from channels chnind

databasef=getenv('DATABASEFILE');
expnum=str2num(getenv('EXPNUM'));
collnum=str2num(getenv('COLLNUM'));

datapath=getenv('DATAPATH');
sdate=getenv('SESSDATE');
samplerate=str2num(getenv('SAMPLERATE'));
trace_b=str2num(getenv('TRACEBLNK')); %trace blank time, ms
% outchans=str2num(getenv('OUTCHANS')); %vector of channels to keep

trwind=str2num(getenv('TRACEWIND')); %window for measuring field
% trjuxt=str2num(getenv('TRACEJUXT')); %juxtalobar amp-lat thresholds
getjuxt=getenv('GETJUXT'); %get juxtalobar line

tracedown=str2num(getenv('TRACEDOWN')); %direction of LFP 1=down, 0=up
tracedlg=getenv('TRACEDLG'); %trace dialog
% trwind=[1 5];
% tracedown=1;
% traceth=str2num(getenv('TRACETH')); %threshold for cleaning noise
tracenum=str2num(getenv('TRACENUM'));    %number of averaged traces
tracevar=getenv('TRACEVAR');   %ordering variable for averaged traces

dtrend=getenv('DETREND');
dtrendt=str2num(getenv('DETRENDT'));
glitchth=str2num(getenv('GLITCHTH'));
glitcht=str2num(getenv('GLITCHT')); %window in trace in which to remove glitches
glitchref=str2num(getenv('GLITCHREF')); %refractory period after glitches

leod=getenv('LEOD');
trjuxt=[];
%channel groups
% changroups=[];
% eval(['changroups=',getenv('CHANGROUPS'),';'])

%spike count window
scwin=str2num(getenv('SCWIN'));

%lpf filter
lp=100*2/samplerate;
[NN, Wn] = buttord( lp, lp .* [0.75], 3, 20);
[B1,A1] = butter(NN,Wn,'high');
t_align=dtrendt-trace_b;
N_align=ceil(t_align*samplerate/1e3);

%% save params to database 
load(databasef);
% experiment(expnum).collection(collnum).segnums=segnums;
experiment(expnum).collection(collnum).tracedown=tracedown;
experiment(expnum).collection(collnum).trwind=trwind;
experiment(expnum).collection(collnum).tracenum=tracenum;
experiment(expnum).collection(collnum).tracevar=tracevar;
experiment(expnum).collection(collnum).dtrend=dtrend;
experiment(expnum).collection(collnum).dtrendt=dtrendt;
experiment(expnum).collection(collnum).glitchth=glitchth;
experiment(expnum).collection(collnum).glitcht=glitcht;
experiment(expnum).collection(collnum).glitchref=glitchref;
experiment(expnum).collection(collnum).leod=leod;
experiment(expnum).collection(collnum).scwin=scwin;
if(strcmp(getjuxt,'off'))
    trjuxt=experiment(expnum).collection(collnum).trjuxt;
elseif(strcmp(getjuxt,'none'))
    trjuxt=nan(16,3);
end
save(databasef,'experiment');
%%
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
amp1=[];
lat1=[];
auc1=[];
sc1=[];
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
    frame.ind0=[frame.ind0;size(frame.data,1)+1];
    withtr=isfield(data.EOD,'traces');
    if(i==1) 
        if(nargin==1)
            file=data.FILE;
            file.offset=[];
            file.title=[];
        end
        if(withtr)
            outchans=[1:size(data.EOD.traces,3)];            
            chnind=outchans;
            for k=1:numel(chnind)
                latnames{k}=['lat',num2str(chnind(k))];
                ampnames{k}=['amp',num2str(chnind(k))];
                aucnames{k}=['auc',num2str(chnind(k))];
            end      
            if(strcmp(tracedlg,'on'))
                FF=figure;
                for k=1:numel(chnind)
                    x=[1:size(data.EOD.traces,2)]/samplerate*1e3+trace_b;
                    y=nanmean(data.EOD.traces(:,:,k));
                    %detrend        
                    if(strcmp(dtrend,'on'))        
                        y=my_detrend(y,N_align);
                    elseif(strcmp(dtrend,'matlab'))        
                        if(numel(N_align))
                            y=detrend(y,'linear',N_align)';
                        else
                            y=detrend(y);
                        end
                    elseif(strcmp(dtrend,'hpf'))        
                        y=filtfilt(B1,A1,y);
                    end                    
                    plot(x,y);
                    [win(k,:),pol(k)]=lfpdlg()
                end
                close(FF);
            else
                if(numel(trwind)>2)
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
            end
            
%             if(numel(trjuxt)==0)
%                trjuxt=zeros(numel(chnind),2);
%             end
        else
           latnames={};
           ampnames={};
           aucnames={};
        end
    end
    if(withtr)
        [lat,amp,auc,avtraces(:,:,:,i)]=get_lfp_stats(data.EOD.traces);
        if(strcmp(getjuxt,'on'))
            load(databasef);
            experiment(expnum).collection(collnum).trjuxt=trjuxt;
            save(databasef,'experiment');
            getjuxt='off';
        end
 
%         gind=get_good_inds(data.EOD.t,data.EOD.traces,amp);        
        gind=1:numel(data.EOD.t);
        lat=lat(gind,:);
        amp=amp(gind,:);
        auc=auc(gind,:);
    else
        avtraces=[];
        lat=[];
        amp=[];
        auc=[];
        gind=1:numel(data.EOD.t);
    end

    deodt=[nan diff(data.EOD.t)];
    eod.t=[eod.t;data.EOD.t(gind)'];
    iei=deodt(gind)';    
    if(isfield(data.EOD,'adcp2p') & strcmp(leod,'on'))
        p=0.1;
        adcp2p=data.EOD.adcp2p(gind,:)/quantile(data.EOD.adcp2p,p);
        adcpr=data.EOD.adcpr(gind,:)/quantile(data.EOD.adcpr,p);
        adcpol=data.EOD.adcpol(gind,:);
        adcnames={'lp2p','lpr','lpol'};
    else
        adcp2p=[];
        adcpr=[];
        adcpol=[];
        adcnames={};
    end
    if(isfield(data.EOD,'raster'))
        u=1;
        for r=1:numel(data.EOD.raster)            
            for p=1:numel(data.EOD.raster{r})                
                for t=1:numel(data.EOD.t)
                    sc(t,u)=sum(data.EOD.raster{r}{p}(:,2)==t & inrange(data.EOD.raster{r}{p}(:,1),scwin));
                end
                scname{u}=['sc',num2str(r),'_',num2str(p)];
                u=u+1;
            end
        end
        sc=sc(gind,:);        
    else
        sc=[];
        scname={};
    end
    lat1=[lat1;lat];
    amp1=[amp1;amp];
    auc1=[auc1;auc];
    sc1=[sc1;sc];
    
    eod.data=[eod.data;iei data.EOD.posture(gind,:)...
        data.EOD.accx(gind)' data.EOD.accy(gind)' data.EOD.accz(gind)' ...
        adcp2p adcpr adcpol];
%          lat amp auc...
%          sc];
    eod.fnames={'iei' data.FILE.model{:} 'accx' 'accy' 'accz' adcnames{:} latnames{:} ampnames{:} aucnames{:} scname{:}};
    eod.ind0=[eod.ind0;size(eod.data,1)+1];
    
    file.offset=[file.offset;data.FILE.offset];    
    
end

eod.data=[eod.data lat1 amp1 auc1 sc1];
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

function [win,pol]=lfpdlg()
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

function [lat1,amp1,auc1,avtraces]=get_lfp_stats(traces)
    x=[1:size(traces,2)]/samplerate*1e3+trace_b;

    for k=1:numel(chnind)
        %clean glitches        
        if(numel(glitchth))        
            t_noise=glitcht;
            N_noise=min(max(ceil(t_noise*samplerate/1e3),1),size(traces,2));            
            ind=find(max(abs(diff(traces(:,N_noise(1):N_noise(2),chnind(k)),1,2)),[],2)>glitchth(k));
            if(numel(ind))
                traces(ind,:,chnind(k))=nan;
            end
        end

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
%         if(~pol(k) & abs(m)<trjuxt(k,1) & (mind/samplerate*1e3+win(k,1))<trjuxt(k,2)) %negative polarity, juxtalobar peak
%             ind=find(x>=trjuxt(k,2) & x<=win(k,2));%narrow window past juxtalobar
%             [m,mind]=min(traces(:,ind,chnind(k)),[],2);
%         end                            
            
        if(pol(k))
%             if(strcmp(getjuxt,'on'))
%                 [trjuxt(k,:)]=get_juxtline(Mind,M);
%                 trjuxt(k,3)=trjuxt(k,3)/samplerate*1e3;%to ms
%             end
%             ijuxt=find(M<trjuxt(k,1)*Mind+trjuxt(k,2)); %juxtalobar peak
%             jind=find(x>=(trjuxt(k,3)+win(k,1)) & x<=win(k,2));%narrow window past juxtalobar
%             [mj,mjind]=max(traces(ijuxt,jind,chnind(k)),[],2);
%             M(ijuxt)=mj;
%             Mind(ijuxt)=mjind+(trjuxt(k,3)*samplerate/1e3);
            amp1(:,k)=M;
%             auc1(:,k)=sum(traces(:,ind,chnind(k)).*(traces(:,ind,chnind(k))>0),2);
            lind=Mind;
        else
            if(strcmp(getjuxt,'on'))
                [trjuxt(k,:)]=get_juxtline(mind,m);
                trjuxt(k,3)=trjuxt(k,3)/samplerate*1e3;%to ms
            end
%             ijuxt=find(abs(m)<trjuxt(k,1) & (mind/samplerate*1e3+win(k,1))<trjuxt(k,2)); %juxtalobar peak
            ijuxt=find(m>trjuxt(k,1)*mind+trjuxt(k,2)); %juxtalobar peak
            jind=find(x>=(trjuxt(k,3)+win(k,1)) & x<=win(k,2));%narrow window past juxtalobar
            [mj,mjind]=min(traces(ijuxt,jind,chnind(k)),[],2);
            m(ijuxt)=mj;
            mind(ijuxt)=mjind+(trjuxt(k,3)-win(k,1))*samplerate/1e3;
            amp1(:,k)=-m;
            auc1(:,k)=-sum(traces(:,ind,chnind(k)).*(traces(:,ind,chnind(k))<0),2);
            lind=mind;            
        end        
        auc1(:,k)=M-m;
%         lind=min(Mind,mind);
        lat1(:,k)=lind/samplerate*1e3+win(k,1);
        I=find(lind==1 | lind==numel(ind)); %no local extremum
        amp1(I,k)=NaN;
        lat1(I,k)=NaN;
        
        if(strcmp(tracevar,'amp'))
            [x,avtraces(:,:,k)]=plot_graded(x,traces(:,:,chnind(k)),[],amp1(:,k),tracenum,0);
        elseif(strcmp(tracevar,'lat'))
            [x,avtraces(:,:,k)]=plot_graded(x,traces(:,:,chnind(k)),[],lat1(:,k),tracenum,0);
        else
            [x,avtraces(:,:,k)]=plot_graded(x,traces(:,:,chnind(k)),[],auc1(:,k),tracenum,0);
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
%     indg=intersect(indg,find(prod(abs(amps)>traceth,2)));
end
   
    function [trj]=get_juxtline(x,y)
        F=figure;
        plot(x,y,'.');
        drawnow;
        [x0,y0]=getline();
        if(numel(x0)==1)
            trj=nan(1,3);
        else
            a=(y0(2)-y0(1))/(x0(2)-x0(1));
            b=y0(1)-a*x0(1);        
            x1=x0(2);
            trj=[a b x1];
        end
        close(F);
    end
end