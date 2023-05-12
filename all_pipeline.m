% set parameters
addpath(genpath('Z:\GitHub\KiloSort')) % path to kilosort folder
addpath(genpath('Z:\GitHub\npy-matlab')) % path to npy-matlab scripts
addpath(genpath('Z:\GitHub\Fish-Pipeline')) % path to pipline
addpath(genpath('Z:\GitHub\Fish-Visualization')) % path to pipline
addpath(genpath('Z:\GitHub\Matlab tools')) % path to pipline

run('.\set_config');
mkdir(ops.spdatpath);
mkdir('metafiles');
output_fname=[ops.analysispath,'\output.mat'];
%% presegment interface
preseg=1;
 if(preseg)
    disp('Pre-segment interface...');
    tic;
    ops=presegment_interface(ops);
    if(exist(output_fname))
        save(output_fname,'ops','-append');
    else
        save(output_fname,'ops');
    end
end

%% analyze segements
segs=1;
if(segs)
    disp('Analyzing Segments....');
    segment_pipeline();        
end

%% spike sort
spsort=1;
if(spsort)
    fs = ops.samplerate; % sampling frequency
    disp('Spike Sorting....');
    for i=1:numel(ops.seg)                
        %  create a channel map file
        Nchannels = numel(cell2mat(ops.seg(i).LFPgroups(ops.seg(i).Spkgroups)));
        if(Nchannels==0)
            continue;
        end
        Nchannels = max(Nchannels,4); %minimum 4 channels for GPU
        connected = true(Nchannels, 1);
        chanMap   = 1:Nchannels;
        chanMap0ind = chanMap - 1;
        kcoords=[];  xcoords=[];    ycoords=[];        
        CHx=[1 ; 2 ; 1 ; 2];
        CHy=[1 ; 1 ; 2 ; 2];
        offset=0;
        for j=1:numel(ops.seg(i).Spkgroups)
            chnum=numel(ops.seg(i).LFPgroups{ops.seg(i).Spkgroups(j)});
            if(numel(ops.seg(i).Spkgroups)==1 & chnum<4)
                chnum=4;
            end
            % grouping of channels (i.e. tetrode groups)
            kcoords   = [kcoords;j*ones(chnum,1)]; 
            xcoords   = [xcoords;  CHx(1:chnum)+offset];
            ycoords   = [ycoords;  CHy(1:chnum)];
            offset=offset+5;
        end        
        save([ops.spdatpath,'\chanMap.mat'], ...
            'chanMap','connected', 'xcoords', 'ycoords', 'kcoords', 'chanMap0ind', 'fs');
        ops.fbinary =[ops.spdatpath,'\spdata_seg',num2str(i),'.bin'];
        if(Nchannels<=8)
            ops.Nfilt= 32; %max(32,4*Nchannels);% number of clusters to use (2-4 times more than Nchan, should be a multiple of 32)     		
        else
            ops.Nfilt= 64;
        end
        %
        tic; % start timer
        %
        if ops.GPU     
            gpuDevice(1); % initialize GPU (will erase any existing GPU arrays)
        end
        [rez, DATA, uproj] = preprocessData(ops); % preprocess data and extract spikes for initialization
        rez= fitTemplates(rez, DATA, uproj);  % fit templates iteratively
        rez= fullMPMU(rez, DATA);% extract final spike times (overlapping extraction)
        % save python results file for Phy
        mkdir(['SortData',num2str(i)]);
        rezToPhy(rez, [ops.analysispath,'\SortData',num2str(i),'\']);
        movefile([ops.spdatpath,'\chanMap.mat'], [ops.analysispath,'\SortData',num2str(i),'\chanMap.mat']);
    end    
end
%% CURATE SORTING USING PHY2
%% precollect interface
precol=1;
load('output.mat', 'ops');
tic;
if(precol)
    disp('Pre-Collection Interface....');
    for i=1:numel(ops.seg) 
        ops=precoll_interface(ops,i);
        save('output.mat','ops','-append');
    end
end
%% after curation: collect data, insert spike data into files
collect_segments();

%% group and PCA
D=load('output.mat','ops','eod');
[ops] = prepca_interface_traces(D.ops,D.eod)
[ops] = indlim_interfact(ops,D.eod)
save('output.mat','ops','-append');
%group and pca
group_pca();    



