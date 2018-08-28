% downscale videos
datapath=getenv('DATAPATH');
sdate=getenv('SESSDATE');
sesspath=[datapath,'\',sdate,'\'];
nfiles=dlmread([sesspath,'counter.txt']); %number of files in recording session
factor=2;
%% loop
for i=0:18
    infile=[sesspath,'video_',num2str(i),'.avi'];
    outfile=[sesspath,'video_',num2str(i),'_scaled.avi'];
    fc=num2str(factor);
    system(['ffmpeg -i ',infile,' -vf scale=iw/',fc,':ih/',fc,' -b:v 1M ',outfile]);
end