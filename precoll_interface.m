function[ops]=precoll_interface(ops,segnum)

%lpf filter
lp=100*2/ops.samplerate;
[NN, Wn] = buttord( lp, lp .* [0.75], 3, 20);
[B1,A1] = butter(NN,Wn,'high');
t_align=ops.detrendt-ops.traceblnk;
N_align=ceil(t_align*ops.samplerate/1e3);

%% get spike data, replace with standard names
% scname={};
% if(numel(ops.seg(segnum).Spkgroups))
%     spk_t=readNPY([ops.analysispath,'\SortData',num2str(segnum),'\spike_times.npy']);
%     spk_c=readNPY([ops.analysispath,'\SortData',num2str(segnum),'\spike_clusters.npy']);
%     cinfo=tdfread([ops.analysispath,'\SortData',num2str(segnum),'\cluster_info.tsv']);    
%     for c=1:numel(cinfo.name)
%         spk_c(spk_c==cinfo.id(c))=cinfo.name(c);
%         scname{c}=['sc',num2str(cinfo.name(c))];
%     end    
%     [spnums,ind]=sort(cinfo.name);
%     scname=scname(ind);
%     ind_start=0;    
% end
%%
    ind=find(ops.seg(segnum).filestat==1);
    filenum=num2str(ops.seg(segnum).files(ind(1)));
    sdate=num2str(ops.seg(segnum).dates(1));
    sesspath=[ops.datapath,'\',sdate,'\'];
    
    s=load([sesspath,'data_',filenum]);
    data=s.data;
    withtr=isfield(data.EOD,'traces');
    chnind=[1:size(data.EOD.traces,3)];            

    get_lfp_stats(data.EOD.traces);
    
    %spike data
%     if(numel(ops.seg(segnum).Spkgroups))
%         ind_end=ind_start+data.FILE.nsamp;
%         spind=find(inrange(spk_t,[ind_start ind_end]));
%         spt=spk_t(spind);
%         spc=spk_c(spind);
%         sc=zeros(numel(data.EOD.t),numel(spnums));
%         for c=1:numel(spnums) %go over all cells
%             sname=spnums(c);
%             rast=double(spt(spc==sname)-ind_start)/ops.samplerate-data.FILE.offset; %spike times in file
%             for e=1:numel(data.EOD.t)   %go over all EODs
%                 t=rast(inrange(rast,[-ops.rasterpre ops.rasterpost]+data.EOD.t(e)))-data.EOD.t(e);
%                 eod.raster{c}=[eod.raster{c};t (eod.ind0(end)+e-1)*ones(size(t))];
%                 sc(e,c)=sum(inrange(t,ops.scwin));
%             end
%         end
% %         eod.raster=[eod.raster;raster];
%     ind_start = ind_end;
%     end
                            


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


function get_lfp_stats(traces)
    x=[1:size(traces,2)]/ops.samplerate*1e3+ops.traceblnk;

    for k=1:numel(chnind)
        %detrend        
        switch ops.detrend
            case 'on'
                traces(:,:,chnind(k))=my_detrend(traces(:,:,chnind(k)),N_align);
            case 'matlab'       
                if(numel(N_align))
                    traces(:,:,chnind(k))=detrend(traces(:,:,chnind(k))','linear',N_align)';
                else
                    traces(:,:,chnind(k))=detrend(traces(:,:,chnind(k))')';
                end
            case 'hpf'
                traces(:,:,chnind(k))=filtfilt(B1,A1,traces(:,:,chnind(k))')';
        end
    end
        
    for i=1:numel(ops.seg(segnum).LFPgroups)
        f(i)=figure;
        for j=1:numel(ops.seg(segnum).LFPgroups{i})
            subplot(2,2,j);
%             H=plot(x,traces(:,:,ops.seg(segnum).LFPgroups{i}(j)));
%             set(H,'Color',0.5*[1 1 1]);
%             hold on;
            H=plot(x,nanmean(traces(:,:,ops.seg(segnum).LFPgroups{i}(j)),1),'k');
            H.LineWidth=3;
            title(num2str(ops.seg(segnum).LFPgroups{i}(j)));
        end
        f(i).Name=['Group ',num2str(i)];
    end
        
    %dialog- choose polarities, windows
    prompt = {'Polarity (1=down):','Windows:'};
    dlgtitle = 'Input';
    dims = [1 35];
    definput = {num2str(ops.seg(segnum).tracedown),num2str(ops.seg(segnum).tracewind)};
    dopts.WindowStyle = 'normal';
    answer = inputdlg(prompt,dlgtitle,dims,definput,dopts);
    close(f(isvalid(f)));

    ops.seg(segnum).tracedown=str2num(answer{1});
    ops.seg(segnum).tracewind=str2num(answer{2});                        
    for j=1:numel(chnind)
        if(numel(ops.seg(segnum).tracewind)>2)
            win(j,:)=ops.seg(segnum).tracewind(j,:);
        else
            win(j,:)=ops.seg(segnum).tracewind;
        end
        if(numel(ops.seg(segnum).tracedown)>1)
            pol(j)=(~ops.seg(segnum).tracedown(j));
        else
            pol(j)=(~ops.seg(segnum).tracedown);
        end
%         if(numel(ops.seg(segnum).getjuxt)>1)
%             gjuxt(j)=ops.seg(segnum).getjuxt(j);
%         else
%             gjuxt(j)=ops.seg(segnum).getjuxt;
%         end
    end
                
    for k=1:numel(chnind)    
        ind=find(x>=win(k,1) & x<=win(k,2));
        [M,Mind]=max(traces(:,ind,chnind(k)),[],2);
        [m,mind]=min(traces(:,ind,chnind(k)),[],2);
            
        if(~pol(k))  %positive 
            [ops.seg(segnum).juxtpoly{k},F]=get_juxttrace(mind,m);
        end
        %remove juxta
        if(numel(ops.seg(segnum).juxtpoly)>=k)
            if(numel(ops.seg(segnum).juxtpoly{k}))
                jind=inpolygon(mind,m,ops.seg(segnum).juxtpoly{k}(:,1),ops.seg(segnum).juxtpoly{k}(:,2));
                if(numel(jind))
                    jtr=nanmean(traces(jind,:,chnind(k))); %compute juxta trace
                    traces(jind,:,chnind(k))=traces(jind,:,chnind(k))-ones(sum(jind),1)*jtr; %remove juxta
                    [m(jind),mind(jind)]=min(traces(jind,ind,chnind(k)),[],2);
                    m(jind)=m(jind)+jtr(mind(jind))';
%                     figure(F);
%                     hold on;
%                     plot(mind,m,'o');
                end
            end
        end
    end
        
end
   
function [pos,F]=get_juxttrace(x,y)
    F=figure;
    plot(x,y,'.');
    drawnow;
    P=drawpolygon();
    pos=P.Position;
    close(F(isvalid(F)));
end
end