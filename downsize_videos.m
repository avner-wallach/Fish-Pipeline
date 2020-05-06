% downscale videos
datapath=getenv('DATAPATH');
sdate=getenv('SESSDATE');
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
%     if(i<=26 | i>=56 | i==45 | i==46)
        system(['ffmpeg -i ',infile,' -vf scale=iw/',fc,':ih/',fc,' -b:v 1M ',outfile]);
%     else
%         system(['ffmpeg -i ',infile,...
%             ' -vf "lutyuv=y=6.47*val-130 [tmp];[tmp] scale=iw/',fc,':ih/',...
%             fc,' [out]" -b:v 1M ',outfile]);
%     end
end