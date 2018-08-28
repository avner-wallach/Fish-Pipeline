function objects_pipline()
%% parameters
% samplerate=2e4;
% 
% M=10; %number of BG images
% K=15;  %median filter kernel size
% online=0; %1: take online tracking data; 0: take offline tracking data
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
BG=[];
objects=[];
corners=[];
for i=1:nfiles
    i
    if(ismember(i-1,skipfiles))
        continue;
    end
    %open video
    vidname=[sesspath,'video_',num2str(i-1),'.avi'];

    if(ismember(i-1,obj_change))
        get_BG(bgframes,vidname);
        get_objects; %plug in hear call to object locating NN                
        save([sesspath,'objects_',num2str(i-1)],'objects','BG','corners');
    end       
end

    function get_objects()      
        F=figure;
        imshow(BG);
        hold on;
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
        [X,Y]=getline(F,'closed');
        if(numel(X)==5)
            corners=[X(1:4) Y(1:4)];
        end
        close(F);
    end    

    function get_BG(imnum,vidname)

        vid=VideoReader(vidname);   
        times=sort(vid.Duration*rand(imnum,1)); %random timepoints
        cdata=zeros(vid.Height,vid.Width,imnum);
        for j=1:imnum
            vid.CurrentTime=times(j);
            F= readFrame(vid);
            FF=mean(F,3);
            cdata(:,:,j)=FF;
        end

        BG=uint8(repmat(median(cdata,3),1,1,3));
    end

        
            
end
