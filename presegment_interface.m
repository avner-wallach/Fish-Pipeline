function ops=presegment_interface(ops)

circle=[54 37 1121 1105];       

for segnum=1:numel(ops.seg)
    
    fi=floor(numel(ops.seg(segnum).files)/2);    
    fprintf('Time %3.0fs. getting file %d of segment %d... \n', toc,fi,segnum);
    filenum=num2str(ops.seg(segnum).files(fi));
    if(numel(ops.seg(segnum).dates)>1)
        sdate=num2str(ops.seg(segnum).dates(fi));
    else
        sdate=num2str(ops.seg(segnum).dates(1));
    end
    
    sesspath=[ops.datapath,'\',sdate,'\'];

    metfilename=[ops.analysispath,'\metafiles\meta_seg',num2str(segnum),'.mat'];
            
    if(~exist(metfilename))
        BG=[];
        objects=[];
        corners=[];
        home=[];
        get_BG;        
        get_objects;
        save([ops.analysispath,'\metafiles\meta_seg',num2str(segnum)],'BG','objects','corners','home','circle');                    
    end
    
    if(numel(ops.chan_num))
        get_amp_data;
    end
end

    function get_amp_data
        [t,amp]=read_bonsai_binary([sesspath,'amp_',filenum],ops.samplerate,ops.chan_num,ops.blocksize,ops.outchans,'adc');
        ind=1:1e5;
        N=numel(ops.seg(segnum).LFPgroups);
        R=ceil(sqrt(N+1));
        C=ceil((N+1)/R);
        F=figure;        
        for n=1:N
            chn=ops.seg(segnum).LFPgroups{n};
            subplot(R,C,n);
            a=diff(amp(:,chn),1,1);
            edges=100:100:3e3;
            for k=1:size(a,2)
                [pks,I]=findpeaks(a(:,k).*(a(:,k)>100),'MinPeakDistance',60);
                h=histc(pks,edges);
                plot(edges,h);
                hold on;
            end
            legend(num2str(chn));
        end
        subplot(R,C,N+1)
        a=max(abs(diff(amp)),[],2);
        edges=100:100:3e3;
        [pks,I]=findpeaks(a.*(a>100),'MinPeakDistance',60);
        h=histc(pks,edges);
        plot(edges,h);
        title(['Max on all channels (EOD chan=0)']);        

        %dialog- choose spk channel
        prompt = {'EOD Chan:','EOD Th:'};
        dlgtitle = 'Input';
        dims = [1 35];
        definput = [{num2str(ops.eodchan),num2str(ops.eodth)}];
        dopts.WindowStyle = 'normal';
        answer = inputdlg(prompt,dlgtitle,dims,definput,dopts);        
        ops.eodchan=str2num(answer{1});
        ops.eodth=str2num(answer{2});
        close(F);     
        
        %get EOD times
        eodind=get_eod_times(t(ind),amp(ind,:),ops);
        
        F=figure;
        plot_raw_traces(amp(ind,:),ops,segnum,eodind,t(ind));
        
%         for n=1:numel(ops.seg(segnum).LFPgroups)
%             F(n)=figure;
%             chn=ops.seg(segnum).LFPgroups{n};
%             subplot(1,2,1);
%             plot(t(ind),amp(ind,chn));            
%             title(['Group:', num2str(n),'  Channels', num2str(chn)]);
%             chname{n}=num2str(chn);
%             pname{n}=['Group ',num2str(n)];
%         end
                
        %dialog- choose channels, spk channel
        prompt ='Spike groups:';
        dlgtitle = 'Input';
        dims = [1 35];
        definput = {num2str(ops.seg(segnum).Spkgroups)};
        dopts.WindowStyle = 'normal';
        answer = inputdlg(prompt,dlgtitle,dims,definput,dopts);
        
        
%         for n=1:numel(ops.seg(segnum).LFPgroups)
%             ops.seg(segnum).LFPgroups{n}=str2num(answer{n});
%         end
        ops.seg(segnum).Spkgroups=str2num(answer{1});

        close(F(isvalid(F)));
        
%         
%         if(~strcmp(ops.eoddetmode,'adc'))
%             
%             if(strcmp(eoddiff,'on'))
%                 a=[0;diff(amp(:,eodchan))];
%             else
%                 a=amp(:,eodchan);
%             end
%         end            
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
                h=imellipse(gca,circle);
                wait(h);
                circle=h.getPosition;
            end                
        end   
        close(F(isvalid(F)));
    end    

            
end