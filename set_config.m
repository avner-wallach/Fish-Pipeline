% setting configuration for experiment pipeline
%% set_config
%database file name
% ops.databasef='Z:\mormyrid_data\analyzed data\database.mat';

% set environemntal variables for pipeline
ops.datapath='Z:\mormyrid_data';
ops.samplerate=30000;
ops.framerate=50;
ops.bits=16;
ops.blocksize=256;
ops.pxlsize=0.423; %pxl to mm conversion

%amp channels
ops.chan_num=16;    %number of amp channels recorded
ops.outchans=[1 1 1 1 1 1 1 0 1 1 1 1 1 1 1 1]; %channels to use

%adc channels
ops.adcchan_num=3; % number of adc channels recorded

%video processing
ops.bgframes=25;    %frames to use for background extraction

ops.framescale=3;   %scaling-down factor for video tracking
ops.trackname='_tracking';  %name extention of tracking file
% ops.trackmode='offline';    %tracking used (on-line/offline)
ops.medfiltk=5;            %tracking median kernel

%eod detection
ops.eoddetmode='amp';   %EOD detection mode amp/adc
ops.eodchan=9;          %channel to detect eod on
ops.eodref=5;           %msec refractory time
ops.eoddiff='on';       %find EOD from amp derivative
ops.eodth=350;          %eod threshold

%filtering
ops.BPF='off';          %bandpass filtering of amp data
ops.hpf_cfr=10;         %HPF cutoff freq
ops.lpf_cfr=3000;       %LPF cutoff freq
ops.n60='off';          %60 hz notch filter

%spike detection
ops.blankgap=1e-3;      %blanking gap
ops.blankpre=1e-3;      %blanking pre-EOD
ops.blankpost=4.5e-3;	%blanking post-EOD
ops.adcthreshold=200;   %ADC noise threshold for blanking 

%raster
ops.rasterpre=0.02;     %s pre EOD
ops.rasterposts=0.04;   %s post EOD
ops.spnum=0;            %running variable for max number of spikes/EOD

%spike count window
ops.scwin=[0.007 0.025];%s post EOD for spike count window

%lfp traces
ops.traceblnk=1;        %trace blank time, ms
ops.tracet=15;          %total trace duration, ms
ops.tracesaved=14;      %saved trace duration

%trace cleaning
ops.traceoffwind=[];%window for removing offsets
ops.detrend='hpf';  %trace de-trending method: matlab=matlab func; on=my func; hpf=high-pass, off=none
ops.detrendt=[];    %for my detrending- timepoints
ops.glitchth=500*ones(size(ops.outchans));%threshold for glitch detection
ops.glicht=[0 5];   %number of consecutive samples under glitch th
ops.glitchref=0;    %refractory period after glitch (150ms in 1st experiment)
ops.traceth=0;      %threshold for cleaning noise

%data collection
ops.tracenum=10;    %number of averaged traces
ops.tracevar='auc'; %ordering variable for averaged traces

%pca coefs
load('Z:\mormyrid_data\analyzed data\pca_coefs.mat');

seg=struct('title',[],'dates',[],'files',[],'blankfiles',[],'blankwins',[],'trackmode','offline','syncfiles',[],...
    'adc_outchans',[0 0 0],'LFPgroups',[],'Spkgroups',[],'tracewind',[],'tracedown',[],'leod',[],'getjuxt',0);
%% segment data

%things valid for all segs
seg.LFPgroups={[1 2 3 4],[5 7],[8 9 10 11],[12 13 14 15]};
seg.Spkgroups=[8 9 10 11];
seg.tracewind=[1 6];
seg.tracedown=1;

ops.seg(1)=seg;
ops.seg(1).title='aqu_brass_leod';
ops.seg(1).dates=20190131;
ops.seg(1).files=[0];
ops.seg(1).blankfiles=[0];
ops.seg(1).blankwins=[3.75 19.5];
ops.seg(1).trackmode='online';
ops.seg(1).adc_outchans=[0 0 1]; 
ops.seg(1).leod=1;

ops.seg(2)=seg;
ops.seg(2).title='aqu_plastic_leod';
ops.seg(2).dates=20190131;
ops.seg(2).files=[1];
ops.seg(2).adc_outchans=[0 0 1]; 
ops.seg(2).leod=1;

ops.seg(3)=ops.seg(2);
ops.seg(3).title='aqu_wood_leod';
ops.seg(3).files=[2];

ops.seg(4)=ops.seg(2);
ops.seg(4).title='aqu_brass';
ops.seg(4).files=[3];

ops.seg(5)=ops.seg(2);
ops.seg(5).title='aqu_plastic';
ops.seg(5).files=[4];

%
ops.seg(6)=seg;
ops.seg(6).title='none';
ops.seg(6).dates=20190131;
ops.seg(6).files=[6:39];

ops.seg(7)=seg;
ops.seg(7).title='none';
ops.seg(7).dates=20190201;
ops.seg(7).files=[0:5];

ops.seg(8)=seg;
ops.seg(8).title='none';
ops.seg(8).dates=20190201;
ops.seg(8).files=[6:10 16:19 25:33];
