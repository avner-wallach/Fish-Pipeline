% downscale videos
datapath=getenv('DATAPATH');
sdate=getenv('SESSDATE');
sesspath=[datapath,'\',sdate,'\'];
nfiles=dlmread([sesspath,'counter.txt']); %number of files in recording session
factor=3;
%% loop
mkdir([sesspath,'\scaled\']);
for i=0:(nfiles-1)
    infile=[sesspath,'video_',num2str(i),'.avi'];
    outfile=[sesspath,'\scaled\video_',num2str(i),'.avi'];
    fc=num2str(factor);
    system(['ffmpeg -i ',infile,' -vf scale=iw/',fc,':ih/',fc,' -b:v 1M ',outfile]);
end