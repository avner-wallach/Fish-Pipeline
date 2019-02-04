function objects_pipline()
%% parameters
datapath=getenv('DATAPATH');
sdate=getenv('SESSDATE');

bgframes=str2num(getenv('BGFRAMES')); %bgframes=100
obj_change=str2num(getenv('OBJ_CHANGE')); %rec. files where object location was changed
skipfiles=str2num(getenv('SKIPFILES')); %files to skip

framescale=str2num(getenv('FRAMESCALE')); %scaling of video file for feature tracking
framecrop=getenv('FRAMECROP');  %cropping of videos for tracking
cropsize=str2num(getenv('CROPSIZE')); %size of cropped image, pre-scaling

sesspath=[datapath,'\',sdate,'\'];
nfiles=dlmread([sesspath,'counter.txt']); %number of files in recording session
%%
for i=1:nfiles
    BG=[];
    objects=[];
    corners=[];
    home=[];
    circle=[];
    
    i
    if(ismember(i-1,skipfiles))
        continue;
    end
    %open video
    vidname=[sesspath,'video_',num2str(i-1),'.avi'];

    if(ismember(i-1,obj_change))
        get_BG(bgframes,vidname);
        get_objects; %plug in hear call to object locating NN                
        save([sesspath,'objects_',num2str(i-1)],'objects','BG','corners','home','circle');
    end       
end

    function get_objects()      
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

    function get_BG(imnum,vidname)

        vid=VideoReader(vidname);   
        times=sort(vid.Duration*rand(imnum,1)); %random timepoints
        cdata=zeros(vid.Height,vid.Width,imnum);
        for j=1:imnum
            vid.CurrentTime=times(j);
exi            FF=mean(F,3);
            cdata(:,:,j)=FF;
        end

        BG=uint8(repmat(median(cdata,3),1,1,3));
    end

        
            
end
