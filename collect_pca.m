function [file,frame,eod] = collect_pca(file,frame,eod)
%COLLECT_PCA get PCA space for posture across all datasets
K=numel(frame)
X=[];
for k=1:K
    X=[X;frame(k).data(:,4:13)];
end
Xm=nanmean(X,1);
[coefs,score,latent,tsq,explained]=pca(X);
for i=1:size(coefs,1)
    pcanames{i}=['pca',num2str(i)];
end

for k=1:K
    %frame data pca
    E=frame(k).data(:,4:13);
    E=E-repmat(Xm,size(E,1),1); %center
    Escore=E*coefs;
    frame(k).data=[frame(k).data Escore];
    frame(k).fnames={frame(k).fnames{:} pcanames{:}};

    %EOD data pca
    E=eod(k).data(:,5:14);
    E=E-repmat(Xm,size(E,1),1); %center
    Escore=E*coefs;
    eod(k).data=[eod(k).data Escore];
    eod(k).fnames={eod(k).fnames{:} pcanames{:}};
    
    file(k).pcacoefs=coefs;
    file(k).explained=explained;
end

end

