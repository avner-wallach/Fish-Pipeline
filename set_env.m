% set environemntal variables for pipeline
setenv('DATAPATH','Z:\AW\mormyrid_data');
setenv('SESSDATE','20180803');
setenv('SAMPLERATE','30000');
setenv('FRAMERATE','50');
setenv('BITS','16');
setenv('BLOCKSIZE','256');
setenv('PXLSIZE','0.423'); %pxl to mm conversion

%amp channels
setenv('CHAN_NUM','4');
setenv('OUTCHANS','1:4')

setenv('BGFRAMES','25');
setenv('SKIPFILES','');
setenv('OBJ_CHANGE','[0]'); %rec. files where object location was changed

setenv('FRAMESCALE','2');
setenv('FRAMECROP','off'); %tracking on cropped image
setenv('CROPSIZE','[800,800]'); %size of cropped image, pre-scaling
% setenv('TRACKNAME','_crpDeepCut_resnet50_fishCroppedScaledAug07shuffle1_1030000');
setenv('TRACKNAME','_scaledDeepCut_resnet50_fishJan30shuffle1_1030000');
setenv('MEDFILTK','5'); %tracking median kernel

%eod detection
setenv('EODCHAN','1');
setenv('EODREF','14.3');%msec refractory time
setenv('EODDIFF','on');%find EOD from amp derivative
setenv('EODTH','500'); %eod threshold

%filtering
setenv('DETREND','off');
setenv('DETRENDT','[0]');
setenv('GLITCHTH','800'); %threshold from glitch detection
setenv('GLITCHT','[0 14]');      %number of consecutive samples under glitch th
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

%trace cleaning
setenv('TRACETH','200'); %threshold for cleaning noise
setenv('TRACEOFFWIND','[]'); %window for removing offsets

%data collection
setenv('TRACEWIND','[2 6;2 6;1.5 4.5;1.5 4.5]'); %window for measuring field
setenv('TRACEDOWN','[0 0 1 1]'); %direction of LFP 1=down, 0=up


