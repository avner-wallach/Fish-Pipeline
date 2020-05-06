%% batch downsize videos
setenv('DATAPATH','Z:\mormyrid_data');
dates={'20200128','20200129'};
for i=1:numel(dates);
    setenv('SESSDATE',dates{i});
    downsize_videos;
end