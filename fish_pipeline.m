function fish_pipeline()
%% parameters
% samplerate=2e4;
% 
% M=10; %number of BG images
% K=15;  %median filter kernel size
% online=0; %1: take online tracking data; 0: take offline tracking data
datapath=getenv('DATAPATH');
sdate=getenv('SESSDATE');
samplerate=str2num(getenv('SAMPLERATE'));
bits=str2num(getenv('BITS'));        %bits=16;
chan_num=str2num(getenv('CHAN_NUM'));%chan_num=16;
blocksize=str2num(getenv('BLOCKSIZE'));%blocksize=256;
eodchan=str2num(getenv('EODCHAN')); %eodchan=3;
bgframes=str2num(getenv('BGFRAMES')); %bgframes=100
obj_change=str2num(getenv('OBJ_CHANGE')); %rec. files where object location was changed
skipfiles=str2num(getenv('SKIPFILES')); %files to skip
framescale=str2num(getenv('FRAMESCALE')); %scaling of video file for feature tracking
sesspath=[datapath,'\',sdate,'\'];
nfiles=dlmread([sesspath,'counter.txt']); %number of files in recording session
%% get adxl parameters
adxl_param=get_adxl_params()
%%
data=[];
for i=1:6%nfiles
    if(ismember(i-1,skipfiles))
        continue;
    end
    %open video
    vidname=[sesspath,'video_',num2str(i-1),'.avi'];
    vid=VideoReader(vidname);
    get_BG(bgframes);
    if(ismember(i-1,obj_change)) %objects change in this video
        get_objects; %plug in hear call to object locating NN
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
        sTS=timestampDecoder(dlmread([sesspath,'video_syncTS',num2str(i-1)]));
        
        pulses=dlmread([sesspath,'pulses_',num2str(i-1)]);
        tsync=pulses(pulses(:,1)==1,2)/samplerate;
        
        %set t=0 to video beginning
        tv0=tv(1);
        tv=tv-tv0;
        sTS=sTS-tv0;    

        %fit sync points
        if(numel(tsync)==(numel(sTS)-1))
            sTS=sTS(1:end-1);            
        elseif(numel(tsync)==(numel(sTS)+1))
            tsync(end)=[];
%         elseif(numel(data.pulses.sync)>(numel(data.sTS)+1))
%             warning(['recording ',num2str(sernum),' sync pulse mismatch']);
%             return;
        end
%         data.samplerate=median(diff(tsync)./diff(sTS));
        
        offset=median(sTS-tsync);
        
        %sync event times
        if(max(pulses(:,1))==2)
            data.EVENT.t=pulses(pulses(:,1)==2,2)/samplerate;
        end
        
        data.FRAME.t=tv;
        data.FILE.offset=offset;
    end

    function get_amp_data(i)
%         R=load([sesspath,'amp_',num2str(i-1)]);
        outchans=[1:16];
        [t,amp]=read_bonsai_data([sesspath,'amp_',num2str(i-1)],samplerate,chan_num,bits,blocksize,outchans);             
        [eodind]=get_eod_times(t,amp);
        [data.EOD.t,data.EOD.lat,data.EOD.amp,data.EOD.traces]=get_lfp_data(t,amp,eodind);
        
    end

    function [eodind]=get_eod_times(t,amp)
        eod_th=0.5; %EOD threshold
        eod_ref=0.01; %EOD refractory time

        ind=find(diff(amp(:,eodchan)>eod_th)==1); %find threshold posedge
        j=2;
        while(j<=numel(ind))
            if((t(ind(j))-t(ind(j-1)))<=eod_ref) %within ref period
                ind(j)=[]; %remove element
            else
                j=j+1;
            end
        end   
        eodind=ind;
%         teod=t(ind);
%         iei=[inf diff(teod)];
    end

    function [teod,lat,amps,traces]=get_lfp_data(t,amp,ind)
        dt=median(diff(t));
        trace_b=0.001; %trace blank time
        trace_T=0.005; %total trace end time
%         trace_t=0.005; %saved trace end time
        trace_B=ceil(trace_b/dt);
        trace_N=ceil(trace_T/dt);
%         trace_n=ceil(trace_t/dt);
        
        ind=ind(ind<(size(amp,1)-trace_N));
        I=ind*ones(1,trace_N) + ones(numel(ind),1)*[0:trace_N-1]; %index matrix
        traces=zeros(numel(ind),trace_N-trace_B,size(amp,2));
        for j=1:size(amp,2)
            a=amp(:,j);
            T=a(I);
            T(:,1:trace_B)=[]; %blank artifact            
            traces(:,:,j)=T;
        end

        [amplitudes,lind]=min(traces,[],2);
        latencies=lind*dt+trace_b;
        % get score- correlation with median trace (to remove 
%         T=median(traces,2);
%         for j=1:size(traces,1)
%             score(i,:)=corr(traces(i,:)',T');
%         end
% 
%         traces=traces(:,1:trace_n);
        amps=permute(amplitudes,[1,3,2]);
        lat=permute(latencies,[1,3,2]);
        teod=t(ind)-data.FILE.offset;
% 
    end            

    function get_aux_data(i)
        fname=[sesspath,'aux_',num2str(i-1)];
        if(exist(fname))
            % accelorometer data
            [aux_t,aux_d]=read_bonsai_data(fname,samplerate/4,3,bits,blocksize/4);             
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
        pitch=atan(y./(sqrt(x.^2 + z.^2)));
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

    function get_BG(imnum)

        times=sort(vid.Duration*rand(imnum,1)); %random timepoints
        cdata=zeros(vid.Height,vid.Width,imnum);
        for j=1:imnum
            j
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
        trackfile=[sesspath,'trackedFeaturesRaw_',num2str(i-1)];
        [num,txt,raw] = xlsread([trackfile,'.csv']);
        [fish,seg,coornames]=get_posture(txt,num);
        data.FRAME.posture=fish(2:end,:)*framescale; %remove 1st frame data
        data.FILE.model=coornames;
        %EOD sample posture
        for j=1:size(fish,2)
            efish(:,j)=interp1(data.FRAME.t,fish(2:end,j)*framescale,data.EOD.t);
        end
        data.EOD.posture=efish;
    end
        
            
end

function adxl_param=get_adxl_params()
if(ismac)
    load('/Users/avner_wallach/Documents/Mormyrid_Data/matlab/ADXL_params_rev1_0.mat');
else
    load('C:\Users\Avner\iCloudDrive\Documents\Mormyrid_Data\matlab\ADXL_params_rev1_0.mat');
end
end