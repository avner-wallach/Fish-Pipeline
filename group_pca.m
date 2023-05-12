function group_pca(segs,infname)
%GROUP_PCA get PCA space for posture across all datasets
%load config
% D=load([pwd,'\metafiles\ops_config.mat']);
% ops=D.ops;
if(nargin<2)
    infname='output.mat';
end
output_fname=['output.mat'];

load(infname);

sdate=num2str(ops.seg(1).dates(1));
archieve_fname=[sdate,'_archieve.mat'];

% save archieved data: file, frame (original + pca'ed), ops
if(isfield(eod,'data'))
    eod=rmfield(eod,'data');
    eod=rmfield(eod,'fnames');
end
save(archieve_fname,'eod','ops');

if(nargin<1)
    segs=1:numel(frame);
end

load([ops.analysispath,'\..\pca_coefs.mat']);

%% posture angles
for i=1:size(tcoefs,1)
    tpcanames{i}=['tpca',num2str(i)];
end

%accelangels
% [acoefs,ascore,alatent,atsq,aexplained]=pca(X(:,8:10));
for i=1:size(acoefs,1)
    apcanames{i}=['apca',num2str(i)];
end
%% move copy data fields
%% frame pca data
for k=segs

    %frame data pca
    if(size(frame(k).data0,2)>5)
        T=frame(k).data0(:,4:10);
        T=T-repmat(Fm(1:7),size(T,1),1); %center
        Tscore=T*tcoefs;
        A=frame(k).data0(:,11:13);
        A=A-repmat(Fm(8:10),size(A,1),1); %center
        Ascore=A*acoefs;

        frame(k).data=[frame(k).data0(:,1:3) Tscore Ascore];
        frame(k).fnames={frame(k).fnames0{1:3} tpcanames{:} apcanames{:}};
    else %headfixed
        frame(k).data=[frame(k).data0(:,1:2)];
        frame(k).fnames={frame(k).fnames0{1:2}};
    end        
end

%% channel pooling names
%%

for k=segs
    k
    cpnames={};
    for j=1:numel(ops.seg(k).LFPgroups);
        cpnames{k}=['lfp',num2str(k)];
    end

    cpnames={[]};
    Xm=[];  Xs=[];
    clear group;
    for i=1:numel(ops.seg(k).LFPgroups_pca)
        i
        cpnames{i}=['lfp',num2str(i)];
        cgroup=ops.seg(k).LFPgroups_pca{i};
        ind=zeros(size(cgroup));
        clear X Y;
        for j=1:numel(cgroup)
            ind(j)=find(cellfun(@(x) strcmp(x,['lfp',num2str(cgroup(j))]),eod(k).fnames0));
            %remove outliers
            I=find(eod(k).CC(:,cgroup(j))>=0);
            I=I(~isoutlier(eod(k).data0(I,ind(j))));
            idx=setdiff(1:size(eod(k).data0,1),I);
            X(:,j)=eod(k).data0(:,ind(j));
            Y(:,j)=eod(k).CC(:,cgroup(j));
            X(idx,j)=nan;
            Y(idx,j)=0;
        end
        [X,xm,xs]=nanzscore(X);
        
        group(:,i)=nansum(X.*Y,2)./nansum(Y,2);        
        Xm=[Xm xm(:)'];
        Xs=[Xs xs(:)'];
%        group(:,i)=nanmedian(medfilt1(nanzscore(eod(k).data0(:,ind)),ops.medfiltk),2);
    end    
    eod(k).LFPm=Xm;
    eod(k).LFPs=Xs;
    
    ind=find(cellfun(@(x) numel(strfind(x,'sc')),eod(k).fnames0));
    sc=eod(k).data0(:,ind);
    scnames=eod(k).fnames0(ind);
    
    %EOD data pca
    if(size(frame(k).data0,2)>5)
        T=eod(k).data0(:,5:11);
        T=T-repmat(Fm(1:7),size(T,1),1); %center
        Tscore=T*tcoefs;
        A=eod(k).data0(:,12:14);
        A=A-repmat(Fm(8:10),size(A,1),1); %center
        Ascore=A*acoefs;    

        eod(k).data=[eod(k).data0(:,1:4) Tscore Ascore group sc];
        eod(k).fnames={eod(k).fnames0{1:4} tpcanames{:} apcanames{:} cpnames{:} scnames{:}};
    else
        eod(k).data=[eod(k).data0(:,1:3) group sc];
        eod(k).fnames={eod(k).fnames0{1:3} cpnames{:} scnames{:}};
    end        
    
    file(k).tcoefs=tcoefs;
    file(k).acoefs=acoefs;
    
end

% save archieved data: file, eod (orginal), frame (original + pca'ed), ops
save(archieve_fname,'file','frame','-append');

% remove original data from eod
eod=rmfield(eod,{'fnames0','data0'});
% save analysis data: file, pca'ed eod, ops
save(output_fname,'file','eod','ops');

end
