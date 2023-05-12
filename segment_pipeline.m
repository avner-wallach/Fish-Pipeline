function segment_pipeline(segs)
%load config
% D=load([pwd,'\metafiles\ops_config.mat']);
tic;

D=load('output.mat','ops');
ops=D.ops;

if(nargin<1)
    segs=1:numel(ops.seg);
end

%% get adxl parameters
adxl_param=get_adxl_params(ops.datapath);
%% filtering coefs
if(strcmp(ops.n60,'on'))
    dnotch = designfilt('bandstopiir','FilterOrder',4, ...
        'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
        'DesignMethod','butter','SampleRate',ops.samplerate);
    [b_dnotch,a_dnotch]=tf(dnotch);
end
% if(strcmp(ops.BPF,'on'))
% dlp = designfilt('bandpassiir','FilterOrder',8, ...
%         'HalfPowerFrequency1',ops.hpf_cfr,'HalfPowerFrequency2',ops.lpf_cfr, ...
%          'SampleRate',ops.samplerate);
%     [b_dlp,a_dlp]=tf(dlp);
% end

%%
data=[];
amp=[];
adc=[];
eodind=[];
for s=1:numel(segs)
    segnum=segs(s);
    if(all(~ops.seg(segnum).filestat)) %unprocessesed segment- remove previous data
        if(exist([ops.spdatpath,'\spdata_seg',num2str(segnum),'.bin']))
            delete([ops.spdatpath,'\spdata_seg',num2str(segnum),'.bin']);
            delete([ops.spdatpath,'\spdata_seg',num2str(segnum),'.log']);
        end
        ops.seg(segnum).eodnum=0;
        ops.seg(segnum).framenum=0;
    end
    for i=1:numel(ops.seg(segnum).files)
        if(ops.seg(segnum).filestat(i)~=0 )  %already collected
            continue;
        end
        
        fprintf('Time %3.0fs. analyzing file %d of segment %d... \n', toc,i,segnum);
        filenum=num2str(ops.seg(segnum).files(i));
        if(numel(ops.seg(segnum).dates)>1)
            sdate=num2str(ops.seg(segnum).dates(i));
        else
            sdate=num2str(ops.seg(segnum).dates);
        end
        sesspath=[ops.datapath,'\',sdate,'\'];

        metfilename=[ops.analysispath,'\metafiles\meta_seg',num2str(segnum),'.mat'];

        if(~exist('BG')) %first file in segment
            if(exist(metfilename))
                load(metfilename);
            elseif(i==1)        
                BG=[];
                objects=[];
                corners=[];
                home=[];
                circle=[];       
                get_BG;        
                get_objects;
                save([ops.analysispath,'\metafiles\meta_seg',num2str(segnum)],'BG','objects','corners','home','circle');                    
            else
                error('run first file first!');
            end
        end

        data.FILE.BG=BG;
        data.FILE.objects=objects;
        data.FILE.corners=corners;
        data.FILE.home=home;
        data.FILE.circle=circle;

        flag=get_sync_data;        
        if(flag==0)
            ops.seg(segnum).filestat(i)=-1; %mark as analyzed
            %save([ops.analysispath,'\metafiles\ops_config.mat'],'ops');
            save([ops.analysispath,'\output.mat'],'ops');
            continue;
        end

        if(numel(ops.outchans))
            get_amp_data;
        end
        if(sum(ops.seg(segnum).adc_outchans)>0)
            get_adc_data;
        end
        if(numel(ops.seg(segnum).Spkgroups))
            save_raw_spfile(amp,ops,segnum,eodind,str2num(filenum));
        end

        get_aux_data;

        if(strcmpi(ops.seg(segnum).trackmode,'offline'))
            get_posture_data;
        else
            get_online_tracking_data;
        end

        save([sesspath,'data_',filenum],'data');
        ops.seg(segnum).filestat(i)=1; %mark as analyzed
        ops.seg(segnum).framenum=ops.seg(segnum).framenum+numel(data.FRAME.t);
        ops.seg(segnum).eodnum=ops.seg(segnum).eodnum+numel(data.EOD.t);        
        %save([ops.analysispath,'\metafiles\ops_config.mat'],'ops');
        save([ops.analysispath,'\output.mat'],'ops','-append');
    end
end

    function flag=get_sync_data
        flag=1;
        v=dlmread([sesspath,'videoTS_',filenum]);
        data.FRAME.counter=cumsum(v(:,2));
        tv=timestampDecoder(v(:,1));

        %video sync timestamps
%         if(ismember(str2num(filenum),ops.seg(segnum).syncfiles))
        if(exist([sesspath,'video_syncTS_offline_',filenum]))
            ind=dlmread([sesspath,'video_syncTS_offline_',filenum]);
            sTS=tv(ind(ind<numel(tv)));
        else
            sTS=timestampDecoder(dlmread([sesspath,'video_syncTS',filenum]));
            if(sTS(1)<tv(1)) %one overflow lost
                sTS=sTS+128;
            end
        end
        
        if(numel(sTS)<3)
            flag=0;
            return;
        end
        
        pulses=dlmread([sesspath,'pulses_',filenum]);
        tsync=pulses(pulses(:,1)==0,2)/ops.samplerate;
        
        %handle video miss-detections
        if(numel(sTS)>numel(tsync)*1.2) %many more video detections
            display('Video sync missdetection, detecting offline...');
            offline_video_sync(ops.datapath,sdate,filenum);
            ind=dlmread([sesspath,'video_syncTS_offline_',filenum]);
            sTS=tv(ind(ind<numel(tv)));
        end            
        
        %set t=0 to video beginning
        tv0=tv(1);
        tv=tv-tv0;
        sTS=sTS-tv0;    

        %pad with nans to get equal lengths
        d=numel(sTS)-numel(tsync);
        if(d>0)            
            tsync=[tsync;nan(d,1)];            
        elseif(d<0)
            sTS=[sTS;nan(-d,1)];
        end
        N=numel(tsync);
        X=[]; Y=[];
        y=tsync;
        for n=1:N
            x=circshift(sTS,n-1);            
            sind=find(abs(diff(x)-diff(y))<1/ops.framerate);
            sind=sind(diff(sind)==1);
%             if(numel(sind)>3 & all(diff(sind)==1))
                X=[X;x(sind+1)]; Y=[Y;y(sind+1)];
%             end
        end
        
        if(numel(X)<3)
            flag=0;
            return;
        end

        c=corr(X,Y);
        if(c>0.75)
            offset=median(X-Y);
        else
            error('small correlation error');
        end

        
%         d=numel(sTS)-numel(tsync);
%         if(d>0)            
%             sTS((end-d+1):end)=[];            
%         elseif(d<0)
%             tsync((end+d+1):end)=[];
%         end
%         %check matching
%         c=corr(diff(sTS),diff(tsync));
%         if(c>0.75)
%             offset=median(sTS-tsync);
%         else
%             X=[]; Y=[];
%             for n=-10:10
%                 x=circshift(sTS,n);
%                 y=tsync;
%                 sind=find(abs(diff(x)-diff(y))<1/ops.framerate);
%                 X=[X;x(sind+1)]; Y=[Y;y(sind+1)];
%             end
%             c=corr(X,Y);
%             if(c>0.75)
%                 offset=median(X-Y);
%             end
% 
%         end
        if(~exist('offset'))
            error('syncpulse error');
        end
        if(abs(offset)>10)
            error('offset error- too large');
        end
        
        %sync event times
        if(max(pulses(:,1))==2)
            data.EVENT.t=pulses(pulses(:,1)==2,2)/ops.samplerate;
        end
        
        data.FRAME.t=tv;
        data.FILE.offset=offset;
    end

    function get_amp_data
        clear t amp;
        if(isfield(ops.seg(segnum),'chan_num') & numel(ops.seg(segnum).chan_num))
            [t,amp]=read_bonsai_binary([sesspath,'amp_',filenum],ops.samplerate,ops.seg(segnum).chan_num,ops.blocksize,ops.seg(segnum).outchans,'adc');
        else
            [t,amp]=read_bonsai_binary([sesspath,'amp_',filenum],ops.samplerate,ops.chan_num,ops.blocksize,ops.outchans,'adc');
        end
        if(ismember(str2num(filenum),ops.seg(segnum).blankfiles))
            amp(inrange(t,ops.seg(segnum).blankwins(ops.seg(segnum).blankfiles==str2num(filenum),:)),:)=0;
        end
        if(isnan(ops.eodth)) %first file
            get_eod_th(amp);
        end        
        amp=filter_voltage(amp);        
        if(~strcmp(ops.eoddetmode,'adc'))
            [eodind]=get_eod_times(t,amp,ops);
            [data.EOD.t,data.EOD.traces]=get_lfp_data(t,amp,eodind);
            data.FILE.nsamp=size(amp,1);
        end            
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
                        
%     function [eodind]=get_eod_times(t,amp)
%         if(ops.eodchan==0)            
%             if(strcmp(ops.eoddiff,'on'))
%                 a=[0;max(abs(diff(amp)),[],2)];
%             else
%                 a=max(amp,[],2);
%             end
%         else
%             if(strcmp(ops.eoddiff,'on'))
%                 a=[0;abs(diff(amp(:,ops.eodchan)))];
%             else
%                 a=amp(:,ops.eodchan);
%             end
%         end            
%         ind=find(diff(a>ops.eodth)==1); %find threshold posedge
%         j=2;
%         while(j<=numel(ind))
%             if((t(ind(j))-t(ind(j-1)))<=ops.eodref/1e3) %within ref period
%                 ind(j)=[]; %remove element
%             else
%                 j=j+1;
%             end
%         end   
%         eodind=ind;
%     end

    function vout=remove_jumps(t,vin,eodinds)        
        for k=1:size(vin,2)
            x=zeros(size(t));
            x(eodinds)=(vin(eodinds+30,k)-vin(eodinds-30,k));
            b=exp(-(1/ops.samplerate)/.12);
            expfilt=filter(1,[1 -b],x);
            vout(:,k)=vin(:,k)-expfilt';
        end
    end

    function [ vin ] = filter_voltage(vin)
%     vout=vin;    
    for j=1:size(vin,2)                        
        %notch
        if(strcmp(ops.n60,'on'))
            tic;
            vin(:,j)=FiltFiltM(b_dnotch,a_dnotch,vin(:,j));
            toc
        end

        if(strcmp(ops.BPF,'on'))
             vin(:,j)=filtfilt(b_dlp,a_dlp,vin(:,j));
        end

    end
    end       

    function [teod,traces]=get_lfp_data(t,amp,ind)
        dt=median(diff(t))*1e3;
        trace_B=ceil(ops.traceblnk/dt);
        trace_N=ceil(ops.tracet/dt);
        trace_M=ceil(ops.tracesaved/dt);
        
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
                
        %save shorter trace
        traces(:,(trace_M+1):end,:)=[];
                
    end            

    function get_aux_data
        fname=[sesspath,'aux_',filenum];
        if(exist(fname))
            % accelorometer data
            [aux_t,aux_d]=read_bonsai_binary(fname,ops.samplerate/4,3,ops.blocksize/4,[1:3],'aux');             
            [t,x,y,z]=get_adxl_vectors(aux_t,aux_d);
            %sample by frame
            data.FRAME.accx=interp1(t,x,data.FRAME.t);
            data.FRAME.accy=interp1(t,y,data.FRAME.t);
            data.FRAME.accz=interp1(t,z,data.FRAME.t);
            %sample by EOD
            if(isfield(data,'EOD'))
                data.EOD.accx=interp1(t,x,data.EOD.t);
                data.EOD.accy=interp1(t,y,data.EOD.t);
                data.EOD.accz=interp1(t,z,data.EOD.t);
            end
        end

    end

    function get_adc_data
        [t,adc]=read_bonsai_binary([sesspath,'adc_',filenum],ops.samplerate,ops.adcchan_num,ops.blocksize,ops.seg(segnum).adc_outchans,'adc');
        if(strcmp(ops.eoddetmode,'adc'))
            [eodind]=get_eod_times(t,adc,ops);
            data.EOD.t=t(eodind)-data.FILE.offset;
            data.FILE.nsamp=size(adc,1);
        end
        [data.EOD.adcp2p,data.EOD.adcpr,data.EOD.adcpol]=get_adc_stats(data.EOD.t,adc);
    end

    function [peak2peak,peakratio,polarity]=get_adc_stats(eodtimes,adc)
        N=2e-3*ops.samplerate; %smaples pre/post eod  
        ind=(eodtimes+data.FILE.offset)*ops.samplerate;        
        I=uint32(ind'*ones(1,2*N+1) + ones(numel(ind),1)*[-N:1:N]); %index matrix 
        I(I<1)=1;
        I(I>size(adc,1))=size(adc,1);
        for m=1:size(adc,2)
%             v=filtfilt(b_dlp,a_dlp,adc(:,m));
            v=adc(:,m);
            traces(:,:,m)=v(I); 
            [mm,im]=min(traces(:,:,m)');
            [M,iM]=max(traces(:,:,m)');
            med=median(traces(:,:,m)');
            peak2peak(:,m)=M-mm;
            peakratio(:,m)=(M-med)./(med-mm);
            polarity(:,m)=(-1).^(iM>im);
        end        
                
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

    function get_posture_data
        %plug in hear call to NN
        get_trackname;
        trackfile=[sesspath,'video_',filenum,ops.trackname];
        [num,txt,raw] = xlsread([trackfile,'.csv']);
        position=[];
        [fish,seg,coornames]=get_posture(txt,num,position,ops.medfiltk);
        
        data.FRAME.posture=fish(2:end,:); %remove 1st frame data
        data.FRAME.posture(:,1:2)=data.FRAME.posture(:,1:2)*ops.framescale;
        data.FILE.model=coornames;
        data.FILE.fish_segment_size=seg*ops.pxlsize*ops.framescale;
        
        %EOD sample posture
        ind=1:min(numel(data.FRAME.t),size(data.FRAME.posture,1));
        for j=1:size(fish,2)
            if(j==3) %for azimuth = unwrap
                efish(:,j)=wrapToPi(interp1(data.FRAME.t(ind),unwrap(data.FRAME.posture(ind,j)),data.EOD.t));
            else
                efish(:,j)=interp1(data.FRAME.t(ind),data.FRAME.posture(ind,j),data.EOD.t);
            end
        end
        data.EOD.posture=efish;
    end

    function get_online_tracking_data
        %plug in here call to NN
        trackfile=[sesspath,'tracking_',filenum];
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

    function get_trackname
        D=dir([sesspath,'*.csv']);
        names={D.name};
        k=0;
        ind=find(cellfun(@(x) numel(strfind(x,['video_0'])),names));       
        while(numel(ind)~=1)            
            k=k+1;
            ind=find(cellfun(@(x) numel(strfind(x,['video_',num2str(k)])),names));       
        end
        ops.trackname=erase(erase(names{ind},['video_',num2str(k)]),'.csv');
    end
        
        
   function get_BG
        vidname=[sesspath,'video_',filenum,'.avi'];
        vid=VideoReader(vidname);   
        times=sort(vid.Duration*rand(ops.bgframes,1)); %random timepoints
        cdata=zeros(vid.Height,vid.Width,ops.bgframes);
        for j=1:ops.bgframes
            vid.CurrentTime=times(j);
            F=readFrame(vid);
            FF=mean(F,3);
            cdata(:,:,j)=FF;
        end
        BG=uint8(repmat(median(cdata,3),1,1,3));
   end

    function get_objects      
        F=figure;
        imshow(BG);
        hold on;
        %get object types
        str = {'Poles','Corners','Home','Circle'};
        [s,v] = listdlg('PromptString','Object Types:',...
                'SelectionMode','multiple',...
                'InitialValue',[1 2],...
                'ListString',str)
        if(v)
            if(ismember(1,s))   %poles
                title('Locate Pole Objects');
                [x,y]=getpts(F);
                ind=find(x>0 & x<size(BG,2) & y>0 & y<size(BG,1));
                for j=1:numel(ind)
                    H=plot(x(ind(j)),y(ind(j)),'x');
                    objects(j).type = questdlg('What kind of object?', ...
                        'Object type', ...
                        'Brass','Plastic','Other','Plastic');
                    objects(j).x=x(ind(j));
                    objects(j).y=y(ind(j));
                    delete(H);
                end
            end
            if(ismember(2,s))   %corners
                title('Locate Corners (rect)');
                h = imrect(gca);
                wait(h);
                rect=h.getPosition;                
                corners(1,:)=rect(1:2);
                corners(2,:)=[rect(1)+rect(3) rect(2)];
                corners(3,:)=[rect(1)+rect(3) rect(2)+rect(4)];
                corners(4,:)=[rect(1) rect(2)+rect(4)];
                delete(h);
            end
            if(ismember(3,s))   %home
                title('Locate Home (points)');
                [X,Y]=getline(F,'closed');
                if(numel(X)==5)
                    home=[X(1:4) Y(1:4)];
                end
            end
            if(ismember(4,s))   %circle
                title('Locate Circle');
                h=imellipse(gca);
                wait(h);
                circle=h.getPosition;
            end                
        end   
        close(F);
    end    

            
end

function adxl_param=get_adxl_params(datapath)
if(ismac)
    load('/Users/avner_wallach/Documents/Mormyrid_Data/matlab/ADXL_params_rev1_0.mat');
else
    load([datapath,'\ADXL_calibrate\ADXL_params_rev1_0.mat']);
end
end
