%% batch downsize videos
setenv('DATAPATH','Z:\mormyrid_data');
dates={'20221116','20221117'};
for i=1:numel(dates);
    setenv('SESSDATE',dates{i});
    downsize_videos;
end