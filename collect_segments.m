function collect_segments(segs)
%load config
D=load('output.mat');
ops=D.ops;

if(nargin<1)
    segs=1:numel(ops.seg);
end

output_fname=[ops.analysispath,'\output.mat'];
tic;

%% load 
if(isfield(D,'frame'))
    %D=load(output_fname);
    frame=D.frame;
    eod=D.eod;
    file=D.file;
end

%lpf filter
lp=100*2/ops.samplerate;
[NN, Wn] = buttord( lp, lp .* [0.75], 3, 20);
[B1,A1] = butter(NN,Wn,'high');
t_align=ops.detrendt-ops.traceblnk;
N_align=ceil(t_align*ops.samplerate/1e3);
%% main loop
for s=segs    
    if(numel(ops.seg(s).Spkgroups) & exist([ops.analysispath,'\SortData',num2str(s),'\cluster_info.tsv']))
        spk_t=readNPY([ops.analysispath,'\SortData',num2str(s),'\spike_times.npy']);
        spk_c=readNPY([ops.analysispath,'\SortData',num2str(s),'\spike_clusters.npy']);
        cinfo=tdfread([ops.analysispath,'\SortData',num2str(s),'\cluster_info.tsv']);    
        if(~isfield(cinfo,'id'))    
            cinfo.id=cinfo.cluster_id;
        end
        spk_n=nan(size(spk_c));
        for c=1:numel(cinfo.name)
            spk_n(spk_c==cinfo.id(c))=cinfo.name(c);
            scname{c}=['sc',num2str(cinfo.name(c))];
        end    
        [spnums,ind]=sort(cinfo.name);
        indu=find(spnums>0);
        spnums=spnums(indu);
        scname=scname(ind(indu));       
        c=numel(spnums);
        %get xcor
        xcor=xcor_all_units(spk_t,spk_n,ops,s); 
        %get ISI histogram
        isi=isi_all_units(spk_t,spk_n,ops,s);
    else
        c=0;
        scname={};
        spnums=[];       
        xcor=[];
    end
    adcp2p=[];
    adcpr=[];
    adcpol=[];
    adcnames={};
    sc=[];            
    
    if(exist('eod') & numel(eod)>=s) %segment already collected
        if(numel(find(ops.seg(s).filestat==2 | ops.seg(s).filestat==-1))==numel(ops.seg(s).filestat)) %all segment analyzed
            continue;                
        end
    end
    if(~exist('eod') | numel(eod)<s | all(ops.seg(s).filestat==1)) %new segement or re-analyze
%         frame(s).t=[];
%         frame(s).data=[];
        frame(s).ind0=[1];
%         frame(s).data0=[];
%         frame(s).fnames={[]};
        frame(s).fnames0={[]};

%         eod(s).t=[];
%         eod(s).data=[];
%         eod(s).data0=[];
        eod(s).ind0=[1];
%         eod(s).fnames={[]};
        eod(s).fnames0={[]};
        eod(s).spk_ind=[0];
%         eod(s).CC=[];

        metfilename=[ops.analysispath,'\metafiles\meta_seg',num2str(s),'.mat'];
        metaf=load(metfilename);        
        file(s).offset=[];
        file(s).title=ops.seg(s).title; 
        file(s).units.xcor=xcor;
        D=fieldnames(metaf);
        for d=1:numel(D)
            file(s).(D{d})=metaf.(D{d});
        end
        
        eod(s).raster=cell(1,c);
        
        lfpnames={};
        lfp=[];
        avtraces=[];
    end        
    
    for i=1:numel(ops.seg(s).files)
        if(ops.seg(s).filestat(i)==2 | ops.seg(s).filestat(i)==-1)  %already collected / not pre-processed
            continue;
        end
        
        fprintf('Time %3.0fs. collecting file %d of segment %d... \n', toc,i,s);
        filenum=num2str(ops.seg(s).files(i));
        if(numel(ops.seg(s).dates)>1)
            sdate=num2str(ops.seg(s).dates(i));
        else
            sdate=num2str(ops.seg(s).dates);
        end
        sesspath=[ops.datapath,'\',sdate,'\'];

        D=load([sesspath,'data_',filenum]);
        data=D.data;
        N=min([numel(data.FRAME.t),size(data.FRAME.posture,1),numel(data.FRAME.accx),numel(data.FRAME.accy),numel(data.FRAME.accz)]);
        framedata=[data.FRAME.posture(1:N,:) data.FRAME.accx(1:N) data.FRAME.accy(1:N) data.FRAME.accz(1:N)];
        if(i==1)
            frame(s).t=nan(ops.seg(s).framenum,1);
            frame(s).data0=nan(ops.seg(s).framenum,size(framedata,2));
            frame(s).counter=0;
            frame(s).fnames0={data.FILE.model{:} 'accx' 'accy' 'accz'};
        end
        frame(s).t(frame(s).counter+(1:N))=data.FRAME.t(1:N);
        frame(s).data0(frame(s).counter+(1:N),:)=framedata;
        frame(s).counter=frame(s).counter+N;        
        frame(s).ind0=[frame(s).ind0;frame(s).counter+1];

        withtr=isfield(data.EOD,'traces');
%         if(i==1) 
%             D=fieldnames(data.FILE);
%             for d=1:numel(D)
%                 file(s).(D{d})=data.FILE.(D{d});
%             end
%         end
        if(withtr)
            chnind=[1:size(data.EOD.traces,3)];                        
            for j=1:numel(chnind)
                lfpnames{j}=['lfp',num2str(chnind(j))];

                if(numel(ops.seg(s).tracewind)>2)
                    win(j,:)=ops.seg(s).tracewind(j,:);
                else
                    win(j,:)=ops.seg(s).tracewind;
                end
                if(numel(ops.seg(s).tracedown)>1)
                    pol(j)=(~ops.seg(s).tracedown(j));
                else
                    pol(j)=(~ops.seg(s).tracedown);
                end
            end
            
            [lfp,avtraces(:,:,:,i),CC]=get_lfp_stats(data.EOD.traces);
        else
           lfpnames={};        
        end

        deodt=[nan diff(data.EOD.t)];
        iei=deodt';    
        if(isfield(data.EOD,'adcp2p') & ops.seg(s).leod)
            p=0.1;
            adcp2p=data.EOD.adcp2p/quantile(data.EOD.adcp2p,p);
            adcpr=data.EOD.adcpr/quantile(data.EOD.adcpr,p);
            adcpol=data.EOD.adcpol;
            adcnames={'lp2p','lpr','lpol'};
        end

        %spike data
        if(c>0)
            ind_start=eod(s).spk_ind(end);
            ind_end=ind_start+data.FILE.nsamp;
            spind=find(inrange(spk_t,[ind_start ind_end]));
            spt=spk_t(spind)-ind_start;
            spc=spk_n(spind);
            sc=zeros(numel(data.EOD.t),numel(spnums));
            for c=1:numel(spnums) %go over all cells
                sname=spnums(c);
%                 rast=double(spt(spc==sname)-ind_start)/ops.samplerate-data.FILE.offset; %spike times in file
                rast=double(spt(spc==sname))/ops.samplerate-data.FILE.offset; %spike times in file;;
                for e=1:numel(data.EOD.t)   %go over all EODs
                    t=rast(inrange(rast,[-ops.rasterpre ops.rasterpost]+data.EOD.t(e)))-data.EOD.t(e);
                    eod(s).raster{c}=[eod(s).raster{c};t (eod(s).ind0(end)+e-1)*ones(size(t))];
                    sc(e,c)=sum(inrange(t,ops.scwin));
                end
            end
            rast=[spt(:) double(spc(:))];
            data.RAST=rast;            
            save([sesspath,'data_',filenum],'data');
            eod(s).spk_ind=[eod(s).spk_ind;ind_end];
    %         eod(s).raster=[eod(s).raster;raster];
            ind_start = ind_end;
        end
        
        eoddata=[iei data.EOD.posture data.EOD.accx' data.EOD.accy' data.EOD.accz' lfp sc adcp2p adcpr adcpol];
        N=size(eoddata,1);
        if(i==1)
            eod(s).t=nan(ops.seg(s).eodnum,1);
            eod(s).data0=nan(ops.seg(s).eodnum,size(eoddata,2));
            eod(s).CC=nan(ops.seg(s).eodnum,size(CC,2));
            eod(s).counter=0;
            eod(s).fnames0={'iei' data.FILE.model{:} 'accx' 'accy' 'accz' lfpnames{:} scname{:} adcnames{:} };
        end        
        eod(s).t(eod(s).counter+(1:N))=data.EOD.t';        
        eod(s).data0(eod(s).counter+(1:N),:)=eoddata;
        eod(s).CC(eod(s).counter+(1:N),:)=CC;
        eod(s).counter=eod(s).counter+N;
        
        eod(s).ind0=[eod(s).ind0;eod(s).counter+1];
        eod(s).avtraces=avtraces;        
        file(s).offset=[file(s).offset;data.FILE.offset];    

        ops.seg(s).filestat(i)=2; %mark as collected
%         save(output_fname,'file','frame','eod','ops');
%        save([ops.analysispath,'\metafiles\ops_config.mat'],'ops');

    end
   save(output_fname,'file','frame','eod','ops');
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

function [lfp,avtraces,CC]=get_lfp_stats(traces)
    x=[1:size(traces,2)]/ops.samplerate*1e3+ops.traceblnk;
    
    for k=1:numel(chnind)
        ind=find(x>=win(k,1) & x<=win(k,2));

        %clean glitches        
        nind=1:size(traces,1);
%         if(numel(ops.glitchth))        
%             N_noise=min(max(ceil(ops.glicht*ops.samplerate/1e3),1),size(traces,2));            
%             wind=find(max(abs(diff(traces(:,ind,chnind(k)),1,2)),[],2)>ops.glitchth(k));
%             if(numel(wind))
%                 traces(wind,:,chnind(k))=nan;
%                 nind=setdiff(nind,wind);
%             end
%         end

        %detrend        
        switch ops.detrend
            case 'on'
                traces(nind,:,chnind(k))=my_detrend(traces(nind,:,chnind(k)),N_align);
            case 'matlab'       
                if(numel(N_align))
                    traces(nind,:,chnind(k))=detrend(traces(nind,:,chnind(k))','linear',N_align)';
                else
                    traces(nind,:,chnind(k))=detrend(traces(nind,:,chnind(k))')';
                end
            case 'hpf'
                traces(nind,:,chnind(k))=filtfilt(B1,A1,traces(nind,:,chnind(k))')';
        end
        
        
        [M,Mind]=max(traces(:,ind,chnind(k)),[],2);
        [m,mind]=min(traces(:,ind,chnind(k)),[],2);
            
        if(pol(k))  %positive 
            amp=M;
            lind=Mind;
        else
%             if(gjuxt(k))
%                 [ops.seg(s).juxtpoly{k}]=get_juxttrace(mind,m);
%                 gjuxt(k)=0;
%             end
            %remove juxta
            if(numel(ops.seg(s).juxtpoly)>=k)
                if(numel(ops.seg(s).juxtpoly{k}))
                    jind=inpolygon(mind,m,ops.seg(s).juxtpoly{k}(:,1),ops.seg(s).juxtpoly{k}(:,2));
                    jtr=nanmean(traces(jind,:,chnind(k)),1); %compute juxta trace
                    traces(jind,:,chnind(k))=traces(jind,:,chnind(k))-ones(sum(jind),1)*jtr; %remove juxta
                    [m(jind),mind(jind)]=min(traces(jind,ind,chnind(k)),[],2);
                    m(jind)=m(jind)+jtr(mind(jind))';
                end
            end
            amp=-m;           
            lind=mind;            
        end        
        p2p=M-m;
        lat=lind/ops.samplerate*1e3+win(k,1);
%         I=find(lind==1 | lind==numel(ind)); %no local extremum
%         amp(I)=NaN;
%         lat(I)=NaN;
%         p2p(I)=NaN;
        switch ops.tracevar
            case 'p2p'
                lfp(:,k)=p2p;
            case 'amp'
                lfp(:,k)=amp;
            case 'lat'
                lfp(:,k)=lat;
        end
        
        %remove outliers
%         indol=find(isoutlier(lfp(:,k)));
        CC(:,k)=corr(traces(:,ind,chnind(k))',nanmedian(traces(:,ind,chnind(k)))');
%         indc=find(CC<-0.25);
%         lfp(indol,k)=nan;
%         lfp(indc,k)=nan;
        
        [x,avtraces(:,:,k)]=plot_graded(x,traces(:,:,chnind(k)),[],lfp(:,k),ops.tracenum,0);        

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
   
    function [pos]=get_juxttrace(x,y)
        F=figure;
        plot(x,y,'.');
        drawnow;
        P=drawpolygon();
        pos=P.Position;
        close(F);
    end
end