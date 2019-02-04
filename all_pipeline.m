% script for running pipeline on all experimental folders
%% set_env
% set environemntal variables for pipeline
setenv('DATAPATH','Z:\mormyrid_data');
setenv('SAMPLERATE','30000');
setenv('FRAMERATE','50');
setenv('BITS','16');
setenv('BLOCKSIZE','256');
setenv('PXLSIZE','0.423'); %pxl to mm conversion

%amp channels
setenv('CHAN_NUM','12');
setenv('OUTCHANS','[1:4 9:12]')

setenv('BGFRAMES','25');

setenv('FRAMESCALE','3');
setenv('FRAMECROP','off'); %tracking on cropped image
setenv('CROPSIZE','[800,800]'); %size of cropped image, pre-scaling
setenv('TRACKNAME','_scaled');
setenv('MEDFILTK','5'); %tracking median kernel

%eod detection
setenv('EODCHAN','1');
setenv('EODREF','14.3');%msec refractory time
setenv('EODDIFF','on');%find EOD from amp derivative
setenv('EODTH','500'); %eod threshold

%filtering
setenv('BPF','off');
setenv('HPF_CFR','10'); %HPF cutoff freq
setenv('LPF_CFR','3000'); %LPF cutoff freq
setenv('N60F','on');  %60 hz notch filter

%lfp traces
setenv('TRACEBLNK','1'); %trace blank time, ms
setenv('TRACET','15'); %total trace duration
setenv('TRACESAVED','14'); %saved trace duration

%video sync
setenv('VIDEO_SYNC','online');
setenv('OFFLINE_SYNC_FILES','[]');

%trace cleaning
setenv('TRACEOFFWIND','[]'); %window for removing offsets
setenv('DETREND','matlab'); %matlab=matlab func; on=my func; off=none
setenv('DETRENDT','[]');
setenv('GLITCHTH','[150 150 150 200 150 300 150 50]'); %threshold from glitch detection
setenv('GLITCHT','[0 5]');      %number of consecutive samples under glitch th
setenv('GLITCHREF','0');      %refractory period after glitch (150ms in 1st experiment)
setenv('TRACETH','0'); %threshold for cleaning noise

%data collection
% setenv('TRACEWIND','[2 6;2 6;1.5 4.5;1.5 4.5]'); %window for measuring field
% setenv('TRACEDOWN','[0 1 1 1]'); %direction of LFP 1=down, 0=up
% setenv('TRACENUM','10');    %number of averaged traces
% setenv('TRACEVAR','amp');   %ordering variable for averaged traces

%%
dates=[20180924 20180925 20180926];
obj_change={'[0 8]','[0 12]','[0 17]'};
channums={'12','12','8'};
outchans={'[1:4 9:12]','[1:4 9:12]','[1:8]'};

skipfiles={'[]','[]','[]'};
offline_sync_files={'[]','[]','[]'};

%% get objects
objects=0;
if(objects)
    for i=1:numel(dates)
        setenv('SESSDATE',num2str(dates(i)));    
        setenv('OBJ_CHANGE',obj_change{i});
        setenv('SKIPFILES',skipfiles{i});
         objects_pipline;
    end
end
%% fish and ephys
fish=0;
if(fish)
    for i=1:numel(dates)
        setenv('SESSDATE',num2str(dates(i)));
        setenv('OBJ_CHANGE',obj_change{i});
        setenv('SKIPFILES',skipfiles{i});
        setenv('CHAN_NUM',channums{i});
        setenv('OUTCHANS',outchans{i})
        setenv('OFFLINE_SYNC_FILES',offline_sync_files{i});
        fish_pipeline;
    end
end
%% collect data
dates={20180924 20180924 20180925 20180925 20180926 20180926};
segs={[0 7],[8 40],[0 11],[12 46],[0 16],[17 27]};
setenv('TRACEWIND','[1 5]'); %window for measuring field
setenv('TRACEDOWN','[1]'); %direction of LFP 1=down, 0=up
setenv('TRACEDLG','off'); %trace dialog
for i=1:numel(dates)
%     frame(i)=struct; eod(i)=struct; file(i)=struct;
    
    for j=1:numel(dates{i})
        setenv('SESSDATE',num2str(dates{i}(j)))
        segnums=[segs{i}(j,1):segs{i}(j,2)];
        if(j==1)
            [frame(i),eod(i),file(i)]=collect_data(segnums);
        else
            [frame(i),eod(i),file(i)]=collect_data(segnums,frame{i},eod{i},file{i});
        end
    end
end

%% pca
% [file,frame,eod] = collect_pca(file,frame,eod);

