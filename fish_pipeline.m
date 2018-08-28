function fish_pipeline()
%% parameters
datapath=getenv('DATAPATH');
sdate=getenv('SESSDATE');
samplerate=str2num(getenv('SAMPLERATE'));
framerate=str2num(getenv('FRAMERATE'));
bits=str2num(getenv('BITS'));        %bits=16;
chan_num=str2num(getenv('CHAN_NUM'));%number of channels recorded 
outchans=str2num(getenv('OUTCHANS')); %vector of channels to keep
blocksize=str2num(getenv('BLOCKSIZE'));%blocksize=256;

dtrend=getenv('DETREND');
dtrendt=str2num(getenv('DETRENDT'));
glitchth=str2num(getenv('GLITCHTH'));
glitcht=str2num(getenv('GLITCHT'));
glitchw=str2num(getenv('GLITCHW'));
bpf=getenv('BPF'); %use bpf
hpf_cfr=str2num(getenv('HPF_CFR')); %HPF cutoff freq
lpf_cfr=str2num(getenv('LPF_CFR')); %LPF cutoff freq
n60f=getenv('N60F');  %60 hz notch filter

eodchan=str2num(getenv('EODCHAN')); %eodchan=3;
eodref=str2num(getenv('EODREF'));%msec refractory time
eoddiff=getenv('EODDIFF');%find EOD from amp derivative
eodth=str2num(getenv('EODTH')); %eod threshold

% bgframes=str2num(getenv('BGFRAMES')); %bgframes=100
obj_change=str2num(getenv('OBJ_CHANGE')); %rec. files where object location was changed
skipfiles=str2num(getenv('SKIPFILES')); %files to skip

framescale=str2num(getenv('FRAMESCALE')); %scaling of video file for feature tracking
framecrop=getenv('FRAMECROP');  %cropping of videos for tracking
cropsize=str2num(getenv('CROPSIZE')); %size of cropped image, pre-scaling

trackname=getenv('TRACKNAME'); %name of tracking file
trace_b=str2num(getenv('TRACEBLNK')); %trace blank time, ms
trace_T=str2num(getenv('TRACET')); %total trace duration
trace_S=str2num(getenv('TRACESAVED')); %saved trace duration

vsync_mode=getenv('VIDEO_SYNC'); 
sesspath=[datapath,'\',sdate,'\'];
nfiles=dlmread([sesspath,'counter.txt']); %number of files in recording session
eod_th=NaN;
%% get adxl parameters
adxl_param=get_adxl_params(datapath);
%%
data=[];
for i=1:nfiles
    i
    if(ismember(i-1,skipfiles))
        continue;
    end
    %open video
    vidname=[sesspath,'video_',num2str(i-1),'.avi'];

    if(ismember(i-1,obj_change) | (~isfield(data,'FILE') & ~exist([sesspath,'data_',num2str(i-2),'.mat'])))%objects change in this video    
        S=load([sesspath,'objects_',num2str(i-1),'.mat']);
        data.FILE.BG=S.BG;
        data.FILE.objects=S.objects;
        data.FILE.corners=S.corners;
%         get_BG(bgframes,vidname);
%         get_objects; %plug in hear call to object locating NN                
    end
    if(~isfield(data,'FILE') & exist([sesspath,'data_',num2str(i-2),'.mat'])) 
        load([sesspath,'data_',num2str(i-2),'.mat']);
    end
%     get_posture(i);

    %     clear data;    
    get_sync_data(i);        
    get_amp_data(i);    
    get_aux_data(i);
    get_posture_data(i);
    
    save([sesspath,'data_',num2str(i-1)],'data');
end


    function get_sync_data(i)
        v=dlmread([sesspath,'videoTS_',num2str(i-1)]);
        data.FRAME.counter=cumsum(v(:,2));
        tv=timestampDecoder(v(:,1));

        %video sync timestamps
        if(strcmp(vsync_mode,'online'))
            sTS=timestampDecoder(dlmread([sesspath,'video_syncTS',num2str(i-1)]));
        else
            ind=dlmread([sesspath,'video_syncTS_offline_',num2str(i-1)]);
            sTS=tv(ind(ind<numel(tv)));
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
            tsync((end-d+1):end)=[];
        end
        %check matching
        c=corr(diff(sTS),diff(tsync))
        if(c>0.75)
            offset=median(sTS-tsync);
        else
            X=[]; Y=[];
            for n=-3:3
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
        
        
        %sync event times
        if(max(pulses(:,1))==2)
            data.EVENT.t=pulses(pulses(:,1)==2,2)/samplerate;
        end
        
        data.FRAME.t=tv;
        data.FILE.offset=offset;
    end

    function get_amp_data(i)
%         outchans=[1:16];
        [t,amp]=read_bonsai_binary([sesspath,'amp_',num2str(i-1)],samplerate,chan_num,blocksize,outchans,'adc');
        if(isnan(eodth)) %first file
            get_eod_th(amp);
        end        
        [eodind]=get_eod_times(t,amp);
        amp=filter_voltage(amp);        
%         amp=remove_jumps(t,amp,eodind);
        [data.EOD.t,data.EOD.traces]=get_lfp_data(t,amp,eodind);
%         [data.EOD.traces,data.EOD.tr_amp,data.EOD.tr_lat]=clean_traces(traces);        
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
            a=[0;diff(amp(:,eodchan))];
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

    function [ vout ] = filter_voltage(vin)
    vout=vin;    
    dnotch = designfilt('bandstopiir','FilterOrder',4, ...
        'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
        'DesignMethod','butter','SampleRate',samplerate);
    if(strcmp(bpf,'on'))
        dlp = designfilt('bandpassiir','FilterOrder',8, ...
                'CutoffFrequency1',hpf_cfr,'CutoffFrequency2',lpf_cfr, ...
                 'SampleRate',samplerate);
    end
    for j=1:size(vin,2)                        
        %notch
        if(strcmp(n60f,'on'))
        vout(:,j)=filtfilt(dnotch,vin(:,j));
        end

        if(strcmp(bpf,'on'))
             vout(:,j)=filtfilt(dlp,vout(:,j));
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
        
        if(strcmp(dtrend,'on'))
            t_align=dtrendt;
            N_align=ceil(t_align/dt);
%             M=mean(traces(:,N_align(1):N_align(2),:),2);
%             traces=traces-repmat(M,1,size(traces,2),1);
            % tr=traces-repmat(median(traces(:,1:N_align,:),2),1,size(traces,2),1);
            for j=1:size(amp,2)
                traces(:,:,j)=detrend(traces(:,:,j)','linear',N_align)';
            end            
        end
        
        if(numel(glitchth))        
            t_noise=glitcht;
            N_noise=min(max(ceil(t_noise/dt),1),size(traces,2));            
            for j=1:size(amp,2)            
                ind=find(max(abs(diff(traces(:,N_noise(1):N_noise(2),j),1,2)),[],2)>glitchth);
                traces(ind,:,j)=nan;
            end
        end
                
        %save shorter trace
        traces(:,(trace_M+1):end,:)=[];
                
    end            

    function get_aux_data(i)
        fname=[sesspath,'aux_',num2str(i-1)];
        if(exist(fname))
            % accelorometer data
            [aux_t,aux_d]=read_bonsai_binary(fname,samplerate/4,3,blocksize/4,[1:3],'aux');             
            [t,roll,pitch,yaw]=get_adxl_angles(aux_t,aux_d);
            %sample by frame
            data.FRAME.roll=interp1(t,roll,data.FRAME.t);
            data.FRAME.pitch=interp1(t,pitch,data.FRAME.t);
            data.FRAME.yaw=interp1(t,yaw,data.FRAME.t);
            %sample by EOD
            data.EOD.roll=interp1(t,roll,data.EOD.t);
            data.EOD.pitch=interp1(t,pitch,data.EOD.t);
            data.EOD.yaw=interp1(t,yaw,data.EOD.t);
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
        if(strcmp(framecrop,'on'))  %video was cropped for tracking
            position=dlmread([sesspath,'tracking_offline_',num2str(i-1)]);
            position(:,1)=position(:,1)-cropsize(1)/2;
            position(:,2)=position(:,2)-cropsize(2)/2;
            position=position/framescale;
            if(size(position,1)<size(num,1))
                num=num(1:end-1,:);
            end
        else
            position=[];
        end
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
        
            
end

function adxl_param=get_adxl_params(datapath)
if(ismac)
    load('/Users/avner_wallach/Documents/Mormyrid_Data/matlab/ADXL_params_rev1_0.mat');
else
    load([datapath,'\ADXL_calibrate\ADXL_params_rev1_0.mat']);
end
end
