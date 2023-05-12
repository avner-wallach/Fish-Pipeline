function fish_pipeline(daynum)
%% parameters

databasef=getenv('DATABASEFILE');
expnum=str2num(getenv('EXPNUM'));

%data
datapath=getenv('DATAPATH');
sdate=getenv('SESSDATE');

%amp parameters
samplerate=str2num(getenv('SAMPLERATE'));
framerate=str2num(getenv('FRAMERATE'));
bits=str2num(getenv('BITS'));        %bits=16;
chan_num=str2num(getenv('CHAN_NUM'));%number of channels recorded 
outchans=str2num(getenv('OUTCHANS')); %vector of channels to keep
blocksize=str2num(getenv('BLOCKSIZE'));%blocksize=256;
blankfile=str2num(getenv('BLANKFILES'));%files in which we blank noise
blankwins=str2num(getenv('BLANKWIN'));%noise blank window

%filtering
bpf=getenv('BPF'); %use bpf
hpf_cfr=str2num(getenv('HPF_CFR')); %HPF cutoff freq
lpf_cfr=str2num(getenv('LPF_CFR')); %LPF cutoff freq
n60f=getenv('N60F');  %60 hz notch filter

%eod detection
eoddetmode=getenv('EODDETMODE'); %adc/amp

eodchan=str2num(getenv('EODCHAN')); %eodchan=3;
eodref=str2num(getenv('EODREF'));%msec refractory time
eoddiff=getenv('EODDIFF');%find EOD from amp derivative
eodth=str2num(getenv('EODTH')); %eod threshold

%channel groups
changroups={};
eval(['changroups=',getenv('CHANGROUPS'),';'])

%spike detection
%blanking
blank_gap=str2num(getenv('BLANKGAP'));
blank_pre=str2num(getenv('BLANKPRE'));
blank_post=str2num(getenv('BLANKPOST'));

%spike detection
reftime=str2num(getenv('REFTIME'));
spikewidth=str2num(getenv('SPIKEWIDTH'));
minpoint=str2num(getenv('MINPOINT'));
thfactor=str2num(getenv('THFACTOR'));
artth=str2num(getenv('ARTTH'));
adcthreshold=str2num(getenv('ADCTHRESHOLD'));
direction=getenv('DIRECTION');

%clustering
blocksize0=str2num(getenv('SBLOCKSIZE')); %min clustering block size
clustnum=str2num(getenv('CLUSTNUM')); %number of clusters in each block
kmeansdims=str2num(getenv('KMEANSDIMS')); %number of clusters in each block

%adc params
adcchan_num=str2num(getenv('ADCCHANNUM')); %total adc channels number
adcchans=str2num(getenv('ADCCHANS')); %adc output channels

% bgframes=str2num(getenv('BGFRAMES')); %bgframes=100
obj_change=str2num(getenv('OBJ_CHANGE')); %rec. files where object location was changed
skipfiles=str2num(getenv('SKIPFILES')); %files to skip

%tracking
framescale=str2num(getenv('FRAMESCALE')); %scaling of video file for feature tracking
% framecrop=getenv('FRAMECROP');  %cropping of videos for tracking
% cropsize=str2num(getenv('CROPSIZE')); %size of cropped image, pre-scaling
tracking_mode=getenv('TRACKMODE'); %ONLINE/OFFLINE
trackname=getenv('TRACKNAME'); %name of tracking file
medfiltk=str2num(getenv('MEDFILTK')); %tracking median kernel

%traces
trace_b=str2num(getenv('TRACEBLNK')); %trace blank time, ms
trace_T=str2num(getenv('TRACET')); %total trace duration
trace_S=str2num(getenv('TRACESAVED')); %saved trace duration

%sync
vsync_mode=getenv('VIDEO_SYNC'); 
offline_sync_files=str2num(getenv('OFFLINE_SYNC_FILES'));

sesspath=[datapath,'\',sdate,'\'];
nfiles=dlmread([sesspath,'counter.txt']); %number of files in recording session
% nfiles=0;
% eod_th=NaN;

dnotch = designfilt('bandstopiir','FilterOrder',4, ...
    'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
    'DesignMethod','butter','SampleRate',samplerate);
[b_dnotch,a_dnotch]=tf(dnotch);
% if(strcmp(bpf,'on'))
    dlp = designfilt('bandpassiir','FilterOrder',8, ...
            'HalfPowerFrequency1',hpf_cfr,'HalfPowerFrequency2',lpf_cfr, ...
             'SampleRate',samplerate);
    [b_dlp,a_dlp]=tf(dlp);

% end

% high-pass for spike extraction
dhpf = designfilt('highpassiir', ...       % Response type
       'StopbandFrequency',200, ...     % Frequency constraints
       'PassbandFrequency',300, ...
       'StopbandAttenuation',55, ...    % Magnitude constraints
       'PassbandRipple',4, ...
       'DesignMethod','cheby1', ...     % Design method
       'MatchExactly','stopband', ...   % Design method options
       'SampleRate',samplerate);               % Sample rate
[b_dhpf,a_dhpf]=tf(dhpf);

%% save params to database
load(databasef);
experiment(expnum).samplerate=samplerate;
experiment(expnum).framerate=framerate;
experiment(expnum).blocksize=blocksize;

experiment(expnum).day(daynum).date=sdate;
experiment(expnum).day(daynum).obj_change=obj_change;
experiment(expnum).day(daynum).chan_num=chan_num;
experiment(expnum).day(daynum).outchans=outchans;
experiment(expnum).day(daynum).blocksize=blocksize;
experiment(expnum).day(daynum).blankfile=blankfile;
experiment(expnum).day(daynum).blankwins=blankwins;
experiment(expnum).day(daynum).bpf=bpf;
experiment(expnum).day(daynum).hpf_cfr=hpf_cfr;
experiment(expnum).day(daynum).lpf_cfr=lpf_cfr;
experiment(expnum).day(daynum).n60f=n60f;
experiment(expnum).day(daynum).eoddetmode=eoddetmode;
experiment(expnum).day(daynum).eodchan=eodchan;
experiment(expnum).day(daynum).eodref=eodref;
experiment(expnum).day(daynum).eoddiff=eoddiff;
experiment(expnum).day(daynum).eodth=eodth;
experiment(expnum).day(daynum).changroups=changroups;
experiment(expnum).day(daynum).blank_gap=blank_gap;
experiment(expnum).day(daynum).blank_pre=blank_pre;
experiment(expnum).day(daynum).blank_post=blank_post;
experiment(expnum).day(daynum).reftime=reftime;
experiment(expnum).day(daynum).spikewidth=spikewidth;
experiment(expnum).day(daynum).minpoint=minpoint;
experiment(expnum).day(daynum).thfactor=thfactor;
experiment(expnum).day(daynum).artth=artth;
experiment(expnum).day(daynum).adcthreshold=adcthreshold;
experiment(expnum).day(daynum).direction=direction;
experiment(expnum).day(daynum).blocksize0=blocksize0;
experiment(expnum).day(daynum).clustnum=clustnum;
experiment(expnum).day(daynum).kmeansdims=kmeansdims;
experiment(expnum).day(daynum).adcchan_num=adcchan_num;
experiment(expnum).day(daynum).adcchans=adcchans;
experiment(expnum).day(daynum).framescale=framescale;
experiment(expnum).day(daynum).tracking_mode=tracking_mode;
experiment(expnum).day(daynum).trackname=trackname;
experiment(expnum).day(daynum).medfiltk=medfiltk;
experiment(expnum).day(daynum).trace_b=trace_b;
experiment(expnum).day(daynum).trace_T=trace_T;
experiment(expnum).day(daynum).trace_S=trace_S;
experiment(expnum).day(daynum).vsync_mode=vsync_mode;
experiment(expnum).day(daynum).offline_sync_files=offline_sync_files;
save(databasef,'experiment');

%% get adxl parameters
adxl_param=get_adxl_params(datapath);
%%
data=[];
amp=[];
adc=[];
for i=1:nfiles
    i
    if(ismember(i-1,skipfiles))
        continue;
    end
    %open video
    vidname=[sesspath,'video_',num2str(i-1),'.avi'];

    if(ismember(i-1,obj_change))% | (~isfield(data,'FILE') & ~exist([sesspath,'data_',num2str(i-2),'.mat'])))%objects change in this video    
        S=load([sesspath,'objects_',num2str(i-1),'.mat']);
        data.FILE.BG=S.BG;
        if(isfield(S,'objects'))
            data.FILE.objects=S.objects;
        else
            data.FILE.objects=[];
        end
        if(isfield(S,'corners'))
            data.FILE.corners=S.corners;
        else
            data.FILE.corners=[];
        end
        if(isfield(S,'home'))
            data.FILE.home=S.home;
        else
            data.FILE.home=[];
        end
        if(isfield(S,'circle'))
            data.FILE.circle=S.circle;
        else
            data.FILE.circle=[];
        end
    elseif(~isfield(data,'FILE'))% & exist([sesspath,'data_',num2str(i-2),'.mat'])) 
        m=2;
        while(~exist([sesspath,'data_',num2str(i-m),'.mat']) & m<i)
            m=m+1;
        end
        if(m==i)
            error('no previous object data!');
        end
        load([sesspath,'data_',num2str(i-m),'.mat']);
    end
    
    get_sync_data(i);        
    
    if(blankfile==i-1)
        blankwin=blankwins(blankfile==(i-1),:);
    else
        blankwin=[nan nan];
    end
    if(numel(outchans))
        get_amp_data(i);
    end
    if(numel(adcchans))
        get_adc_data(i);
    end
    if(numel(changroups))
        [data.SPIKES]=get_spike_times(amp,changroups,eodind,b_dhpf,a_dhpf,data.FILE.offset,adc,sesspath,daynum);
    end
    
    get_aux_data(i);

    if(strcmpi(tracking_mode,'OFFLINE'))
        get_posture_data(i);
    else
        get_online_tracking_data(i);
    end
        
    save([sesspath,'data_',num2str(i-1)],'data');
end


    function get_sync_data(i)
        v=dlmread([sesspath,'videoTS_',num2str(i-1)]);
        data.FRAME.counter=cumsum(v(:,2));
        tv=timestampDecoder(v(:,1));

        %video sync timestamps
        if(strcmp(vsync_mode,'offline') | ismember(i-1,offline_sync_files))
            ind=dlmread([sesspath,'video_syncTS_offline_',num2str(i-1)]);
            sTS=tv(ind(ind<numel(tv)));
        else
            sTS=timestampDecoder(dlmread([sesspath,'video_syncTS',num2str(i-1)]));
            if(sTS(1)<tv(1)) %one overflow lost
                sTS=sTS+128;
            end
        end
        
        pulses=dlmread([sesspath,'pulses_',num2str(i-1)]);
        tsync=pulses(pulses(:,1)==1,2)/samplerate;
        
        %set t=0 to video beginning
        tv0=tv(1);
        tv=tv-tv0;
        sTS=sTS-tv0;    

        %fit sync points
        d=numel(sTS)-numel(tsync);
        if(d>0)            
            sTS((end-d+1):end)=[];            
        elseif(d<0)
            tsync((end+d+1):end)=[];
        end
        %check matching
        c=corr(diff(sTS),diff(tsync))
        if(c>0.75)
            offset=median(sTS-tsync);
        else
            X=[]; Y=[];
            for n=-10:10
                x=circshift(sTS,n);
                y=tsync;
                sind=find(abs(diff(x)-diff(y))<1/framerate);
                X=[X;x(sind+1)]; Y=[Y;y(sind+1)];
            end
            c=corr(X,Y);
            if(c>0.75)
                offset=median(X-Y);
            end

        end
        if(~exist('offset'))
            error('syncpulse error');
        end
        if(abs(offset)>10)
            error('offset error- too large');
        end

        
        %sync event times
        if(max(pulses(:,1))==2)
            data.EVENT.t=pulses(pulses(:,1)==2,2)/samplerate;
        end
        
        data.FRAME.t=tv;
        data.FILE.offset=offset;
    end

    function get_amp_data(i)
        [t,amp]=read_bonsai_binary([sesspath,'amp_',num2str(i-1)],samplerate,chan_num,blocksize,outchans,'adc');
        amp(inrange(t,blankwin),:)=0;
        if(isnan(eodth)) %first file
            get_eod_th(amp);
        end        
        amp=filter_voltage(amp);        
        if(~strcmp(eoddetmode,'adc'))
            [eodind]=get_eod_times(t,amp);
            [data.EOD.t,data.EOD.traces]=get_lfp_data(t,amp,eodind);
        end
%         generate_amp_file(amp,changroups,eodind,b_dhpf,a_dhpf); %save filtered data in one big file for FAST
%         end
            
    end

    function get_eod_th(amp)
        if(strcmp(eoddiff,'on'))
            a=[0;diff(amp(:,eodchan))];
        else
            a=amp(:,eodchan);
        end
        
        F=figure;
        A=axis;
        ind=1:1e5;
        H=plot(a(ind));
        prompt = {'EOD Threshold:'};
        dlg_title = '';
        num_lines = 1;
        defaultans = {'5000'};
        opt.WindowStyle='normal';
        answer = inputdlg(prompt,dlg_title,num_lines,defaultans,opt);
        eodth=str2num(answer{1});
        close(F);
    end
                        
    function [eodind]=get_eod_times(t,amp)
        if(strcmp(eoddiff,'on'))
            a=[0;abs(diff(amp(:,eodchan)))];
        else
            a=amp(:,eodchan);
        end
        ind=find(diff(a>eodth)==1); %find threshold posedge
        j=2;
        while(j<=numel(ind))
            if((t(ind(j))-t(ind(j-1)))<=eodref/1e3) %within ref period
                ind(j)=[]; %remove element
            else
                j=j+1;
            end
        end   
        eodind=ind;
    end

    function vout=remove_jumps(t,vin,eodinds)        
        for k=1:size(vin,2)
            x=zeros(size(t));
            x(eodinds)=(vin(eodinds+30,k)-vin(eodinds-30,k));
            b=exp(-(1/samplerate)/.12);
            expfilt=filter(1,[1 -b],x);
            vout(:,k)=vin(:,k)-expfilt';
        end
    end

    function [ vin ] = filter_voltage(vin)
%     vout=vin;    
    for j=1:size(vin,2)                        
        %notch
        if(strcmp(n60f,'on'))
            tic;
            vin(:,j)=FiltFiltM(b_dnotch,a_dnotch,vin(:,j));
            toc
        end

        if(strcmp(bpf,'on'))
             vin(:,j)=filtfilt(b_dlp,a_dlp,vin(:,j));
        end

    end
    end       

    function [teod,traces]=get_lfp_data(t,amp,ind)
        dt=median(diff(t))*1e3;
        trace_B=ceil(trace_b/dt);
        trace_N=ceil(trace_T/dt);
        trace_M=ceil(trace_S/dt);
        
        ind=ind(ind<(size(amp,1)-trace_N));
        I=ind*ones(1,trace_N) + ones(numel(ind),1)*[0:trace_N-1]; %index matrix
        traces=zeros(numel(ind),trace_N-trace_B,size(amp,2));
        for j=1:size(amp,2)
            a=amp(:,j);
            T=a(I);
            T(:,1:trace_B)=[]; %blank artifact            
            traces(:,:,j)=T;
        end
        
        teod=t(ind)-data.FILE.offset;
        
%         if(strcmp(dtrend,'on'))
%             t_align=dtrendt;
%             N_align=ceil(t_align/dt);
% %             M=mean(traces(:,N_align(1):N_align(2),:),2);
% %             traces=traces-repmat(M,1,size(traces,2),1);
%             % tr=traces-repmat(median(traces(:,1:N_align,:),2),1,size(traces,2),1);
%             for j=1:size(amp,2)
%                 traces(:,:,j)=detrend(traces(:,:,j)','linear',N_align)';
%             end            
%         end
        
%         if(numel(glitchth))        
%             t_noise=glitcht;
%             N_noise=min(max(ceil(t_noise/dt),1),size(traces,2));            
%             for j=1:size(amp,2)            
%                 ind=find(max(abs(diff(traces(:,N_noise(1):N_noise(2),j),1,2)),[],2)>glitchth);
%                 traces(ind,:,j)=nan;
%             end
%         end
                
        %save shorter trace
        traces(:,(trace_M+1):end,:)=[];
                
    end            

    function get_aux_data(i)
        fname=[sesspath,'aux_',num2str(i-1)];
        if(exist(fname))
            % accelorometer data
            [aux_t,aux_d]=read_bonsai_binary(fname,samplerate/4,3,blocksize/4,[1:3],'aux');             
%             [t,roll,pitch,yaw]=get_adxl_angles(aux_t,aux_d);
            [t,x,y,z]=get_adxl_vectors(aux_t,aux_d);
            %sample by frame
%             data.FRAME.roll=interp1(t,roll,data.FRAME.t);
%             data.FRAME.pitch=interp1(t,pitch,data.FRAME.t);
%             data.FRAME.yaw=interp1(t,yaw,data.FRAME.t);
            data.FRAME.accx=interp1(t,x,data.FRAME.t);
            data.FRAME.accy=interp1(t,y,data.FRAME.t);
            data.FRAME.accz=interp1(t,z,data.FRAME.t);
            %sample by EOD
%             data.EOD.roll=interp1(t,roll,data.EOD.t);
%             data.EOD.pitch=interp1(t,pitch,data.EOD.t);
%             data.EOD.yaw=interp1(t,yaw,data.EOD.t);
            if(isfield(data,'EOD'))
                data.EOD.accx=interp1(t,x,data.EOD.t);
                data.EOD.accy=interp1(t,y,data.EOD.t);
                data.EOD.accz=interp1(t,z,data.EOD.t);
            end
        end

    end

    function get_adc_data(i)
        [t,adc]=read_bonsai_binary([sesspath,'adc_',num2str(i-1)],samplerate,adcchan_num,blocksize,adcchans,'adc');
        if(strcmp(eoddetmode,'adc'))
            [eodind]=get_eod_times(t,adc);
            data.EOD.t=t(eodind)-data.FILE.offset
        end
        [data.EOD.adcp2p,data.EOD.adcpr,data.EOD.adcpol]=get_adc_stats(data.EOD.t,t,adc);
    end

    function [peak2peak,peakratio,polarity]=get_adc_stats(eodtimes,t,adc)
        N=2e-3*samplerate; %smaples pre/post eod  
        ind=(eodtimes+data.FILE.offset)*samplerate;        
        I=uint32(ind'*ones(1,2*N+1) + ones(numel(ind),1)*[-N:1:N]); %index matrix 
        I(I<1)=1;
        I(I>size(adc,1))=size(adc,1);
        for m=1:size(adc,2)
            v=filtfilt(b_dlp,a_dlp,adc(:,m));
            traces(:,:,m)=v(I); 
            [mm,im]=min(traces(:,:,m)');
            [M,iM]=max(traces(:,:,m)');
            med=median(traces(:,:,m)');
            peak2peak(:,m)=M-mm;
            peakratio(:,m)=(M-med)./(med-mm);
            polarity(:,m)=(-1).^(iM>im);
        end        
                
    end            

    function [t,roll,pitch,yaw]=get_adxl_angles(aux_t,aux_d)

        D=50; %decimate factor

        %convert to g
        x=(aux_d(:,1)-adxl_param.x0)/adxl_param.xg; %left=positive
        y=(aux_d(:,2)-adxl_param.y0)/adxl_param.yg; %front=positive
        z=(aux_d(:,3)-adxl_param.z0)/adxl_param.zg; %down=positive

        %total acceleration
        A=sqrt(x.^2 + y.^2 + z.^2);
        roll=atan(x./(sqrt(y.^2 + z.^2)));
        pitch=atan(y./(sqrt(x.^2 + z.^2)))-pi/2;
        yaw=atan(z./(sqrt(x.^2 + y.^2)));

        %down sample
        roll=decimate(roll,D);
        pitch=decimate(pitch,D);
        yaw=decimate(yaw,D);
        t=decimate(aux_t-data.FILE.offset,D);

    end

    function [t,x,y,z]=get_adxl_vectors(aux_t,aux_d)

        D=50; %decimate factor

        %convert to g
        x=(aux_d(:,1)-adxl_param.x0)/adxl_param.xg; %left=positive
        y=(aux_d(:,2)-adxl_param.y0)/adxl_param.yg; %front=positive
        z=(aux_d(:,3)-adxl_param.z0)/adxl_param.zg; %down=positive

        %down sample
        x=decimate(x,D);
        y=decimate(y,D);
        z=decimate(z,D);
        t=decimate(aux_t-data.FILE.offset,D);

    end

    function get_objects()      
        F=figure;
        imshow(data.FILE.BG);
        hold on;
        [x,y]=getpts(F);
        ind=find(x>0 & x<size(data.FILE.BG,2) & y>0 & y<size(data.FILE.BG,1));
        for j=1:numel(ind)
            H=plot(x(ind(j)),y(ind(j)),'x');
            data.FILE.objects(j).type = questdlg('What kind of object?', ...
                'Object type', ...
                'Brass','Plastic','Other','Plastic');
            data.FILE.objects(j).x=x(ind(j));
            data.FILE.objects(j).y=y(ind(j));
            delete(H);
        end
        close(F);
    end    

    function get_BG(imnum,vidname)

        vid=VideoReader(vidname);   
        times=sort(vid.Duration*rand(imnum,1)); %random timepoints
        cdata=zeros(vid.Height,vid.Width,imnum);
        for j=1:imnum
%             j
            vid.CurrentTime=times(j);
            F= readFrame(vid);
            FF=mean(F,3);
        %     FF=mean(F,3)';
        %     Q=1082*918;
        %     FF=reshape(FF(1:Q),918,1082)';
            cdata(:,:,j)=FF;
        end

        data.FILE.BG=uint8(repmat(median(cdata,3),1,1,3));
    end

    function get_posture_data(i)
        %plug in hear call to NN
        trackfile=[sesspath,'video_',num2str(i-1),trackname];
        [num,txt,raw] = xlsread([trackfile,'.csv']);
%         if(strcmp(framecrop,'on'))  %video was cropped for tracking
%             position=dlmread([sesspath,'tracking_offline_',num2str(i-1)]);
%             position(:,1)=position(:,1)-cropsize(1)/2;
%             position(:,2)=position(:,2)-cropsize(2)/2;
%             position=position/framescale;
%             if(size(position,1)<size(num,1))
%                 num=num(1:end-1,:);
%             end
%         else
            position=[];
%         end
        [fish,seg,coornames]=get_posture(txt,num,position);
        
        %interpolate isolated nans
%         for k=1:size(position,2)
%             ind=find(isnan(position(2:end-1,k)))+1;
%             position(ind,k)=(position(ind-1,k)+position(ind+1,k))/2;
%         end
        
        data.FRAME.posture=fish(2:end,:); %remove 1st frame data
        data.FRAME.posture(:,1:2)=data.FRAME.posture(:,1:2)*framescale;
        data.FILE.model=coornames;
        %EOD sample posture
        ind=1:min(numel(data.FRAME.t),size(data.FRAME.posture,1));
        for j=1:size(fish,2)
            efish(:,j)=interp1(data.FRAME.t(ind),data.FRAME.posture(ind,j),data.EOD.t);
        end
        data.EOD.posture=efish;
    end

    function get_online_tracking_data(i)
        %plug in hear call to NN
        trackfile=[sesspath,'tracking_',num2str(i-1)];
        a_x_y=dlmread(trackfile);
        data.FRAME.posture=a_x_y(:,2:3);
        
        %EOD sample posture
        ind=1:min(numel(data.FRAME.t),size(data.FRAME.posture,1));
        for j=1:size(data.FRAME.posture,2)
            efish(:,j)=interp1(data.FRAME.t(ind),data.FRAME.posture(ind,j),data.EOD.t);
        end
        data.EOD.posture=efish;
        data.FILE.model={'X','Y'};
    end

            
end

function adxl_param=get_adxl_params(datapath)
if(ismac)
    load('/Users/avner_wallach/Documents/Mormyrid_Data/matlab/ADXL_params_rev1_0.mat');
else
    load([datapath,'\ADXL_calibrate\ADXL_params_rev1_0.mat']);
end
end
