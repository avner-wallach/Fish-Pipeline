function[ops]=pregroup_interface(segs)
%load config
D=load([pwd,'\metafiles\ops_config.mat']);
ops=D.ops;

if(nargin<1)
    segs=1:numel(ops.seg);
end

%load data
output_fname=[ops.analysispath,'\output.mat'];
load(output_fname);


%% 

for i=1:numel(segs)
    X=eod(i).data0;
    for j=1:numel(ops.seg(i).LFPgroups)
        clear cind;
        for k=1:numel(ops.seg(i).LFPgroups{j})
            cind(k)=find(cellfun(@(x) strcmp(x,['lfp',num2str(ops.seg(i).LFPgroups{j}(k))]),eod(i).fnames0));
        end
        x=nanzscore(X(:,cind));
        F(j)=plot_group(x,eod(i).fnames0(cind));
        pname{j}=['Group ',num2str(j)];
        chname{j}=num2str(ops.seg(i).LFPgroups{j});
    end
    %dialog- choose channels, spk channel
    prompt = [pname];
    dlgtitle = 'Input';
    dims = [1 35];
    definput = [chname,{''}];
    dopts.WindowStyle = 'normal';
    answer = inputdlg(prompt,dlgtitle,dims,definput,dopts);

    for n=1:numel(ops.seg(i).LFPgroups)
        ops.seg(i).LFPgroups{n}=str2num(answer{n});
    end
    close(F(isvalid(F)));
end

end

function F=plot_group(y,names)
F=figure;
N=size(y,2);
nm={};
for n=1:(N-1)
    for m=(n+1):N
        plot(y(:,n),y(:,m),'.');
        nm=[nm [names{n},'-',names{m}]];
        hold on;
    end
end
legend(nm,'Location','northwest');
set(gca,'Xlim',[-5 5],'Ylim',[-5 5]);

%inset
for n=1:N
    ind=setdiff(1:N,n);
    S(n)=nanmean(nanstd(y(:,ind),[],2));
end
A=axes;
A.Position=[0.6661    0.1600    0.2354    0.2590];
bar(1:N,S);
A.XTickLabel=names;

end
        