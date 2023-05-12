function [ops] = prepca_interface(ops,eod)
%PREPCA_INTERFACE go over all segments, check for failed LFP and revise
%   LFPgroups

for j=1:numel(ops.seg(1).LFPgroups)
    for i=1:numel(ops.seg)
        lfgorup=ops.seg(i).LFPgroups{j};
        cols=[];
        for k=1:numel(lfgorup)
            cols=[cols find(cellfun(@(x) (strcmp(x,['lfp',num2str(lfgorup(k))])),eod(i).fnames0))];
            if(isfield(eod(i),'CC'))
                ind{k}=find(eod(i).CC(:,lfgorup(k))>=0);
                ind{k}=ind{k}(find(~isoutlier(eod(i).data0(ind{k},cols(k)))));
            else
                ind{k}=1:size(eod(i).data0,1);
            end
        end
        F(i)=figure;
        if(numel(cols)==2)
            inds=intersect(ind{1},ind{2});
            plot(nanzscore(eod(i).data0(inds,cols(1))),nanzscore(eod(i).data0(inds,cols(2))),'.');
        else
            for k=1:numel(cols)
                plot(nanzscore(eod(i).data0(ind{k},cols(k))));
                hold on;            
            end
        end
        lfpname{i}=num2str(lfgorup);
    end            
    ops=lfpdlg(lfpname,ops);
    close(F(:)); 
    
    %get ind limits for spike rate
    K=150;
    for i=1:numel(ops.seg)
        cols=find(cellfun(@(x) numel(strfind(x,'sc')),eod(i).fnames0));
        C=numel(cols);
        F=figure;
        F.Name=['segment ',num2str(i)];
        M=ceil(sqrt(C));
        N=ceil(C/M);
        for k=1:C
            subplot(N,M,k);
            plot(smooth(eod(i).data0(:,cols(k)),K));
            title(eod(i).fnames0(cols(k)));
        end
        ops=sratedlg(ops,eod(i).fnames(cols),i);
        close(F);
    end

end

function [ops]=lfpdlg(chans,ops)
    prompt = cell(numel(chans),1);
    dlg_title = 'LFP channels to take';
    num_lines = 1;
    defaultans = chans;
    dopts.WindowStyle = 'normal';    
    answers = inputdlg(prompt,dlg_title,num_lines,defaultans,dopts);
    for i=1:numel(chans)
        ops.seg(i).LFPgroups_pca{j}=str2num(answers{i});
    end
end

function ops=sratedlg(ops,scnames,segnum)
    prompt = cell(numel(scnames),1);
    dlg_title = 'Spike Rate Limiters';
    num_lines = numel(scnames);
    defaultans = repmat([1 inf],num_lines,1);
    dopts.WindowStyle = 'normal';    
    answers = inputdlg(prompt,dlg_title,num_lines,defaultans,dopts);
    for n=1:num_lines
        ops.seg(segnum).ind_lim(n,:)=str2num(answers{n});
    end
end

end
