% downscale videos
datapath=getenv('DATAPATH');
sdate=getenv('SESSDATE');
% datapath = 'Z:\mormyrid_data';
% sdate = '20210309';
sesspath=[datapath,'\',sdate,'\'];
nfiles=dlmread([sesspath,'counter.txt']); %number of files in recording session
% nfiles=2;
factor=3;
fc=num2str(factor);
%% loop
mkdir([sesspath,'\scaled\']);
for i=0:(nfiles-1)
    infile=[sesspath,'video_',num2str(i),'.avi'];
    outfile=[sesspath,'\scaled\video_',num2str(i),'.avi'];
    system(['ffmpeg -i ',infile,' -vf scale=iw/',fc,':ih/',fc,' -b:v 1M ',outfile]);
end