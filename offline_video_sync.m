function offline_video_sync(datapath,datename,filenum)
% detect and log sync pulses from video file
% ROI=[1045 1156 1075 1200];
% ROI=[1150 1127 1157 1135];
ROI = [1042 1155 1088 1186];
roith=250;
actth=10;
ind=[];

vidname=[datapath,filesep,datename,filesep,'video_',filenum,'.avi']
vid=VideoReader(vidname);
k=1;
while(vid.hasFrame)
    F=vid.readFrame;
%     if(k==1)
%         figure;
%         imshow(F);
%         h=drawrectangle;
%         ROI=round(h.Position);
%         ROI(3)=h.Position(1)+h.Position(3);
%         ROI(4)=h.Position(2)+h.Position(4);
%     end
    roi=F(ROI(2):ROI(4),ROI(1):ROI(3),1);
    rois(k)=max(roi(:));
    roim(k)=mean(roi(:));
    act(k)=sum(roi(:)>=roith);
    k=k+1;
end
% ind=find(diff(act>actth)==1);
ind=find(diff(roim>40)==1);
if(numel(ind)<100)
    ind
end
fname=[datapath,filesep,datename,filesep,'video_syncTS_offline_',filenum];
dlmwrite(fname,ind);
end