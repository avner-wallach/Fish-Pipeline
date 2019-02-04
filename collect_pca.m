function [pframe,peod,file] = collect_pca(file,frame,eod,chanidx,chanpool)
%COLLECT_PCA get PCA space for posture across all datasets
K=numel(frame)

%% posture angles
%concat all frame posture 
X=[];
for k=1:K
    X=[X;frame(k).data(:,4:13)];
end
Xm=nanmean(X,1);
%tracking angels
[tcoefs,tscore,tlatent,ttsq,texplained]=pca(X(:,1:7));
for i=1:size(tcoefs,1)
    tpcanames{i}=['tpca',num2str(i)];
end

%accelangels
[acoefs,ascore,alatent,atsq,aexplained]=pca(X(:,8:10));
for i=1:size(acoefs,1)
    apcanames{i}=['apca',num2str(i)];
end

%% frame pca data
for k=1:K
    %frame data pca
    T=frame(k).data(:,4:10);
    T=T-repmat(Xm(1:7),size(T,1),1); %center
    Tscore=T*tcoefs;
    A=frame(k).data(:,11:13);
    A=A-repmat(Xm(8:10),size(A,1),1); %center
    Ascore=A*acoefs;
        
    pframe(k).t=frame(k).t;
    pframe(k).ind0=frame(k).ind0;
    pframe(k).data=[frame(k).data(:,1:3) Tscore Ascore];
    pframe(k).fnames={frame(k).fnames{1:3} tpcanames{:} apcanames{:}};
end

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
    %EOD data pca
    T=eod(k).data(:,5:11);
    T=T-repmat(Xm(1:7),size(T,1),1); %center
    Tscore=T*tcoefs;
    A=eod(k).data(:,12:14);
    A=A-repmat(Xm(8:10),size(A,1),1); %center
    Ascore=A*acoefs;    
    
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
    peod(k).data=[eod(k).data(:,1:4) Tscore Ascore E pool];
    peod(k).fnames={eod(k).fnames{1:4} tpcanames{:} apcanames{:} eod(k).fnames{chanidx} cpnames{:}};
    
    file(k).tcoefs=tcoefs;
    file(k).texplained=texplained;
    file(k).acoefs=acoefs;
    file(k).aexplained=aexplained;
    
end

end

