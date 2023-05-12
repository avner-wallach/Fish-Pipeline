function [ops] = indlim_interfact(ops,eod)
%   get ind limits for spike rate

K=150;
for i=1:numel(ops.seg)
    if(isfield(eod(i),'fnames'))
        fnames=eod(i).fnames;
        data=eod(i).data;
    else
        fnames=eod(i).fnames0;
        data=eod(i).data0;
    end
    cols=find(cellfun(@(x) numel(strfind(x,'sc')),fnames));
    C=numel(cols);
    if(~C)
        continue;
    end
    F=figure;
    F.Name=['segment ',num2str(i)];
    M=ceil(sqrt(C));
    N=ceil(C/M);
    for k=1:C
        subplot(N,M,k);
        plot(smooth(data(:,cols(k)),K));
        title(fnames(cols(k)));
    end
    ops=sratedlg(ops,fnames(cols),i);
    close(F);
end

function ops=sratedlg(ops,scnames,segnum)
    prompt = cell(numel(scnames),1);
    dlg_title = 'Spike Rate Limiters';
    num_lines =1;% numel(scnames);
    defaultans = cell(numel(scnames),1);
    for i=1:numel(scnames)
        if(isfield(ops.seg(segnum),'ind_lim') & numel(ops.seg(segnum).ind_lim>0))
            defaultans{i}=num2str(ops.seg(segnum).ind_lim(i,:));
        else
            defaultans{i}='1 inf';
        end
    end
    dopts.WindowStyle = 'normal';    
    answers = inputdlg(prompt,dlg_title,num_lines,defaultans,dopts);
    for n=1:numel(scnames)
        ops.seg(segnum).ind_lim(n,:)=str2num(answers{n});
    end
end

end

