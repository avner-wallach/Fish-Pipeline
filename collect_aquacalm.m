function [peod] = collect_aquacalm(file,eod,chanidx,chanpool)
%COLLECT_PCA get PCA space for posture across all datasets
K=numel(eod)

%% eod data
X=[];
for k=1:K
    X=[X;eod(k).data(:,chanidx)];
end
Ym=nanmean(X,1);
Ys=nanstd(X);
for i=1:numel(chanidx)
    if(strfind(eod(1).fnames{chanidx(i)},'lat'))
        Yd(i)=-1;
    else
        Yd(i)=1;
    end
end
%% channel pooling names
cpnames={};
for k=1:numel(chanpool)
    cpnames{k}=['pool',num2str(k)];
end
%%

for k=1:K

    E=eod(k).data(:,chanidx);
    Ym1=nanmean(E,1);
    Ys1=nanstd(E);
    E=E-repmat(Ym1,size(E,1),1); %center
    E=E./repmat(Ys1,size(E,1),1).*repmat(Yd,size(E,1),1); %flip latency
    
    %pool
    clear pool;
    pool=[];
    for i=1:numel(chanpool)
        pool(:,i)=nanmedian(E(:,chanpool{i}),2);
    end

    peod(k).t=eod(k).t;    
    peod(k).ind0=eod(k).ind0;    
    peod(k).data=[eod(k).data(:,1:4) pool];
    peod(k).fnames={eod(k).fnames{1:4}  cpnames{:}};
    
    
end

end

