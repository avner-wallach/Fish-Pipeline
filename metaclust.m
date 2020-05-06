function [ ] = metaclust(daynum)
%METACLUST meta clustering of centroids
databasef=getenv('DATABASEFILE');
expnum=str2num(getenv('EXPNUM'));

%% parameters
datapath=getenv('DATAPATH');
sdate=getenv('SESSDATE');
sesspath=[datapath,'\',sdate,'\'];
nfiles=dlmread([sesspath,'counter.txt']); %number of files in recording session
skipfiles=str2num(getenv('SKIPFILES')); %files to skip
rasterpre=str2num(getenv('RASTERTPRE'));
rasterpost=str2num(getenv('RASTERTPOST'));

%% save to database
load(databasef);
experiment(expnum).day(daynum).clust_skipfiles=skipfiles;
experiment(expnum).day(daynum).rasterpre=rasterpre;
experiment(expnum).day(daynum).rasterpost=rasterpost
save(databasef,'experiment');

%% get all centroid and isihist data
if(~exist([sesspath,'clustind.mat']))
    for i=1:nfiles
        i
        if(ismember(i-1,skipfiles))
            continue;
        end

        data=getfield(load([sesspath,'data_',num2str(i-1)]),'data');

        G=numel(data.SPIKES);    %number of electrode groups
        if(~exist('centroids'))    
            centroids=cell(1,G);
            isihist=cell(1,G);
            S=zeros(1,G);
        end        
        for g=1:G
            B=numel(data.SPIKES(g).blocks);
            for b=1:B
               centroids{g}=[centroids{g};data.SPIKES(g).blocks(b).centroids];
               isihist{g}=[isihist{g};data.SPIKES(g).blocks(b).isihist];
            end
            numblocks(i,g)=B;        
        end
    end
    isibins=data.SPIKES(1).isibins;
    for g=1:G        
        indices{g}=ones(size(centroids{g},1),1);    
        % perform meta clustering
        cluststr=MyPCA(centroids{g},isihist{g},isibins,indices{g});
        indices{g}=cluststr.indices;
        save([sesspath,'clustind'],'indices');
        S(g)=max(S(g),max(indices{g}));
    end
else
    load([sesspath,'clustind.mat']);
    for g=1:numel(indices)
        S(g)=max(indices{g});
    end
end
%% associate spikes with meta-clusters
T=10; %ms max lag for xcor
k=ones(size(S));
for i=1:nfiles
    i
    if(ismember(i-1,skipfiles))
        continue;
    end

    data=getfield(load([sesspath,'data_',num2str(i-1)]),'data');
    
    G=numel(data.SPIKES);    %number of electrode groups
    for g=1:G
        cents=cell(S(g),1); %collecting centroids per unit        
        spind=[];
        B=numel(data.SPIKES(g).blocks);
        for b=1:B
            lookup=indices{g}(k(g):k(g)+size(data.SPIKES(g).blocks(b).centroids,1)-1); %metaclusters associated with centroids in this block
            for s=1:S(g)
                cents{s}=[cents{s};data.SPIKES(g).blocks(b).centroids(lookup==s,:)];
            end
            spind=[spind;lookup(data.SPIKES(g).blocks(b).indices)]; %associate spikes with metaclusters
            k(g)=k(g)+size(data.SPIKES(g).blocks(b).centroids,1);
        end
        data.SPIKES(g).raster=[data.SPIKES(g).raster(:,1) spind];        

        for s1=1:S(g) %for each unit
            %compute centroid stats
            data.SPIKES(g).trace(s1).m=nanmean(cents{s1},1);
            data.SPIKES(g).trace(s1).s=nanstd(cents{s1},[],1);
            %compute xcor
            for s2=s1:S
                [c,bins]=MyXcor(data.SPIKES(g).raster(data.SPIKES(g).raster(:,2)==s1,1)*1e3,...
                    data.SPIKES(g).raster(data.SPIKES(g).raster(:,2)==s2,1)*1e3,T);
                data.SPIKES(g).xcor(:,s1,s2)=c;
            end
            data.SPIKES(g).xbins=bins;
        end
    end
    
    [data.EOD.raster] = get_eod_rasters( data,S );
    
    save([sesspath,'data_',num2str(i-1)],'data');
end

