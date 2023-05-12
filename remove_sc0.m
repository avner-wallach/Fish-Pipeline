%remove 0 unit
fname='20190624_archieve.mat';
load(fname);
for s=1:numel(eod)
    ind0=find(cellfun(@(x) strcmp(x,'sc0'),eod(s).fnames0));
    eod(s).fnames0(ind0)=[];
    eod(s).data0(:,ind0)=[];
    eod(s).raster(1)=[];
end
save(fname,'eod','file','frame','ops');

load('output.mat');
for s=1:numel(eod)
    ind0=find(cellfun(@(x) strcmp(x,'sc0'),eod(s).fnames));
    eod(s).fnames(ind0)=[];
    eod(s).data(:,ind0)=[];
    eod(s).raster(1)=[];
end
save('output.mat','eod','file','ops');
    