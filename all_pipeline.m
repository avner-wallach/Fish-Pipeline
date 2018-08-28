% script for running pipeline on all experimental folders
set_env;
dates=[20180803 20180804 20180806];
segs=[0 18;0 48;0 10];
% segs=[0 1;0 1;1 2];
trackname={'_crpDeepCut_resnet50_fishCroppedScaledAug07shuffle1_1030000',...
    '_scaledDeepCut_resnet50_fishJan30shuffle1_1030000',...
    '_scaledDeepCut_resnet50_fishJan30shuffle1_1030000'};
obj_change={'0','[0 46]','[0]'};
cropping={'on','off','off'};

%get objects
for i=1:numel(dates)
    setenv('SESSDATE',num2str(dates(i)));    
    setenv('OBJ_CHANGE',obj_change{i});
%     objects_pipline;
end

%fish and ephys
for i=1:numel(dates)
    setenv('FRAMECROP',cropping{i});
    setenv('SESSDATE',num2str(dates(i)));
    setenv('TRACKNAME',trackname{i});
    setenv('OBJ_CHANGE',obj_change{i});
%     fish_pipeline;
end

%collect data
for i=1:numel(dates)
    setenv('SESSDATE',num2str(dates(i)));
    segnums=[segs(i,1):segs(i,2)];
    [frame(i),eod(i),file(i)]=collect_data(segnums);
end

%pca
[file,frame,eod] = collect_pca(file,frame,eod);

