function [pframe,peod,file] = collect_pca(file,frame,eod,groupfeat,tcoefs,acoefs,Fm)
%COLLECT_PCA get PCA space for posture across all datasets
eval(['changroups=',getenv('CHANGROUPS'),';'])
% eval(['groupfeat=',getenv('GROUPFEAT'),';']);
% tcoefs=str2num(getenv('TCOEFS'));
% acoefs=str2num(getenv('ACOEFS'));
databasef=getenv('DATABASEFILE');
expnum=str2num(getenv('EXPNUM'));
collnum=str2num(getenv('COLLNUM'));

K=numel(frame)
%% save params to database
load(databasef);
experiment(expnum).collection(collnum).changroups=changroups;
experiment(expnum).collection(collnum).groupfeat=groupfeat;
save(databasef,'experiment');

%% posture angles
%concat all frame posture 
% X=[];
% for k=1:K
%     X=[X;frame(k).data(:,4:13)];
% end
% Xm=nanmean(X,1);
%tracking angels
% [tcoefs,tscore,tlatent,ttsq,texplained]=pca(X(:,1:7));
for i=1:size(tcoefs,1)
    tpcanames{i}=['tpca',num2str(i)];
end

%accelangels
% [acoefs,ascore,alatent,atsq,aexplained]=pca(X(:,8:10));
for i=1:size(acoefs,1)
    apcanames{i}=['apca',num2str(i)];
end

%% frame pca data
for k=1:K
    %frame data pca
    if(size(frame(k).data,2)>5)
        T=frame(k).data(:,4:10);
        T=T-repmat(Fm(1:7),size(T,1),1); %center
        Tscore=T*tcoefs;
        A=frame(k).data(:,11:13);
        A=A-repmat(Fm(8:10),size(A,1),1); %center
        Ascore=A*acoefs;

        pframe(k).t=frame(k).t;
        pframe(k).ind0=frame(k).ind0;
        pframe(k).data=[frame(k).data(:,1:3) Tscore Ascore];
        pframe(k).fnames={frame(k).fnames{1:3} tpcanames{:} apcanames{:}};
    else %headfixed
        pframe(k).t=frame(k).t;
        pframe(k).ind0=frame(k).ind0;
        pframe(k).data=[frame(k).data(:,1:2)];
        pframe(k).fnames={frame(k).fnames{1:2}};
    end
        
end

%% eod data

% X=[];
% for k=1:K
%     X=[X;eod(k).data(:,chanidx)];
% end
% Ym=nanmean(X,1);
% Ys=nanstd(X);
% for i=1:numel(chanidx)
%     if(strfind(eod(1).fnames{chanidx(i)},'lat'))
%         Yd(i)=-1;
%     else
%         Yd(i)=1;
%     end
% end
%% channel pooling names
cpnames={};
for k=1:numel(changroups)
    cpnames{k}=['group',num2str(k)];
end
%%

for k=1:K
    for i=1:numel(changroups)
        cgroup=changroups{i};
        ind=zeros(size(cgroup));
        for j=1:numel(cgroup)
            ind(j)=find(cellfun(@(x) strcmp(x,[groupfeat{i},num2str(cgroup(j))]),eod(k).fnames));
        end
        group(:,i)=nanmedian(nanzscore(eod(k).data(:,ind)),2);
    end    
    ind=find(cellfun(@(x) numel(strfind(x,'sc')),eod(k).fnames));
    sc=eod(k).data(:,ind);
    scnames=eod(k).fnames(ind);
    
    %EOD data pca
    if(size(frame(k).data,2)>5)
        T=eod(k).data(:,5:11);
        T=T-repmat(Fm(1:7),size(T,1),1); %center
        Tscore=T*tcoefs;
        A=eod(k).data(:,12:14);
        A=A-repmat(Fm(8:10),size(A,1),1); %center
        Ascore=A*acoefs;    

        peod(k).t=eod(k).t;    
        peod(k).ind0=eod(k).ind0;    
        peod(k).data=[eod(k).data(:,1:4) Tscore Ascore group sc];
        peod(k).fnames={eod(k).fnames{1:4} tpcanames{:} apcanames{:} cpnames{:} scnames{:}};
    else
        peod(k).t=eod(k).t;    
        peod(k).ind0=eod(k).ind0;    
        peod(k).data=[eod(k).data(:,1:3) group sc];
        peod(k).fnames={eod(k).fnames{1:3} cpnames{:} scnames{:}};
    end        
    
    file(k).tcoefs=tcoefs;
%     file(k).texplained=texplained;
    file(k).acoefs=acoefs;
%     file(k).aexplained=aexplained;
    file(k).changroups=changroups;
    
end

end

