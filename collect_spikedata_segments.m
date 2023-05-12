function collect_spikedata_segments(segs)
%load config
load('output.mat','file','ops');

if(nargin<1)
    segs=1:numel(ops.seg);
end

output_fname=[ops.analysispath,'\output.mat'];
tic;

%% load 
%% main loop
for s=segs    
    if(numel(ops.seg(s).Spkgroups) & exist([ops.analysispath,'\SortData',num2str(s),'\cluster_info.tsv']))
        spk_t=readNPY([ops.analysispath,'\SortData',num2str(s),'\spike_times.npy']);
        spk_c=readNPY([ops.analysispath,'\SortData',num2str(s),'\spike_clusters.npy']);
        cinfo=tdfread([ops.analysispath,'\SortData',num2str(s),'\cluster_info.tsv']);    
        if(~isfield(cinfo,'id'))    
            cinfo.id=cinfo.cluster_id;
        end
        spk_n=nan(size(spk_c));
        for c=1:numel(cinfo.name)
            spk_n(spk_c==cinfo.id(c))=cinfo.name(c);
            scname{c}=['sc',num2str(cinfo.name(c))];
        end    
        [spnums,ind]=sort(cinfo.name);
        indu=find(spnums>0);
        spnums=spnums(indu);
        scname=scname(ind(indu));       
        c=numel(spnums);
        %get xcor
        %xcor=xcor_all_units(spk_t,spk_n,ops,s); 
        %get ISI histogram
        isi=isi_all_units(spk_t,spk_n,ops,s);
    else
%         xcor=[];
        isi=[];
    end
%     file(s).units.xcor=xcor;
    file(s).units.isi=isi;
   save(output_fname,'file','-append');

end