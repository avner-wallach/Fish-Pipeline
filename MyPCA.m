function output_struct=MyPCA(centroids,isihist,isibins,indices)

% spike sorter version for efish data
clear data clusters;
global clusters handles params data outstr flag;
data.traces=centroids;     %segmented threshold crossings
data.isihist=isihist;
data.isibins=isibins;
clusters=indices;
flag=0;
%% define structures
% parameters struct
params.samplerate=30; %khz
params.feats=[1 2];
% params.spikewidth=2;
% params.windowwidth=.5;
% params.windowstart=0.9;
% params.minpoint=1;
% params.reftime=1;
% params.thfactor=4.5;
% params.direction='both';
% params.art_reftime=.8;
% params.art_width=1;
% params.art_thresh=1e4;
% params.art_smoothorder=20;
% params.art_N=10;   
% params.art_Fc=1.25; %1Khz highpass filtering
% params.art_direction='rising';
% params.blankback=13;
% params.blanklength=25;
% params.lambda=0.9; %minimal posterior for classification

%analysis params

% data struct
data.score=[];      %PCA score per event
data.coeff=[];      %PCA coefficients

%% design GUI
handles.figure1 = figure;

% Create axes
handles.axes1 = axes('Parent',handles.figure1,...
    'Position',[0.025 0.12 0.575 0.8]);
box(handles.axes1,'on');
hold(handles.axes1,'all');

% Create title
title({'PCA'});

% Create axes
handles.axes2=axes('Parent',handles.figure1,...
    'Position',[0.64 0.55 0.33 0.36]);

% Create title
title({'Traces'});

% Create axes
handles.axes3=axes('Parent',handles.figure1,...
    'Position',[0.64 0.12 0.33 0.36]);

% Create title
title({'ISI Histograms'});

x=10;
y=20;
w=150;
% handle1 = uicontrol('Style','pushbutton','String','Choose Files',...
%     'Position',[x y w 50],'Callback',{@choose_files});
% x=x+w;
handle2 = uicontrol('Style','pushbutton','String','Set Parameters',...
    'Position',[x y w 50],'Callback',{@set_params});
x=x+w;
% handle2 = uicontrol('Style','pushbutton','String','Set PCA Window',...
%     'Position',[x y w 50],'Callback',{@set_pcawind});
% x=x+w;
handle7 = uicontrol('Style','pushbutton','String','Auto Cluster',...
    'Position',[x y w 50],'Callback',{@auto_clust});
x=x+w;
handle5 = uicontrol('Style','pushbutton','String','Merge Clusters',...
    'Position',[x y w 50],'Callback',{@merge_clusters});
x=x+w;
handle7 = uicontrol('Style','pushbutton','String','Split Cluster',...
    'Position',[x y w 50],'Callback',{@auto_split});
x=x+w;
handle4 = uicontrol('Style','pushbutton','String','Split Outliers',...
    'Position',[x y w 50],'Callback',{@outlier_split});
x=x+w;
% handle6 = uicontrol('Style','pushbutton','String','Plot Corr',...
%     'Position',[x y w 50],'Callback',{@draw_xcor});
% x=x+w;
% handle6 = uicontrol('Style','pushbutton','String','Show Trace',...
%     'Position',[x y w 50],'Callback',{@show_trace});
% x=x+w;
handle6 = uicontrol('Style','pushbutton','String','Show Shapes',...
    'Position',[x y w 50],'Callback',{@show_shapes});
x=x+w;
handle8 = uicontrol('Style','pushbutton','String','Remove Artifact',...
    'Position',[x y w 50],'Callback',{@remove_artifact});
x=x+w;
handle100 = uicontrol('Style','pushbutton','String','END',...
     'Position',[x y w 50],'Callback',{@export_data});

load_data();

while(~flag)
    pause(1);
end
output_struct=outstr;

end

%% event functions
function set_params(hObj,event)
    global handles params data;

    prompt={'Features'};
    name='Update Parameters';
    numlines=1;
    options.Resize='on';
    options.WindowStyle='normal';
    defaultanswer={''};
    answer1=inputdlg(prompt,name,numlines,defaultanswer,options);
    if(numel(answer1))
        if(numel(str2num(answer1{1})))
            params.feats=str2num(answer1{1});
        end

    end
    
    draw_clusters;

end
    
function load_data()
    global clusters handles params data;
%     data.times=[];      %event times
    data.score=[];      %PCA score per event
    data.coeff=[];      %PCA coefficients
    data.art_traces=[]; %artifact traces
%     data.art_times=[];  %artifact times

%     clusters=ones(size(data.traces,1),1);

    disp('Computing PCA.....');
    [data.coeff,data.score,latent]=pca(data.traces);
    data.eigvals=cumsum(latent)./sum(latent);
    
    data.t=(1:size(data.traces,2))/params.samplerate;

    axes(handles.axes1);
    cla;

    draw_clusters;

end

function auto_clust(hObj,event)
    global clusters handles data params;    

    answer=inputdlg({'Number of Clusters?'});
    if(numel(answer) & numel(str2num(answer{1})))
        clust_num=str2num(answer{1});
        display('automatic clustering...');
        clust_idx=kmeans(data.score(clusters~=0,params.feats),clust_num); %clust_idx is without artifact
        clusters(clusters~=0)=clust_idx;
        
        draw_clusters;
%         draw_xcor;
    end

end

function auto_split(hObj,event)
    global clusters handles data params;    

    answer=inputdlg({'Which Cluster to split?','Number of Sub-Clusters?','features?'});
    if(numel(answer) & numel(str2num(answer{1})))
        split_clust=str2num(answer{1});
        clust_num=str2num(answer{2});
        feats=str2num(answer{3});
        display('automatic sub-clustering...');
        idx=find(clusters==split_clust);
        if(numel(idx))
            if(numel(feats))
                %compute PCA to cluster
                [coeff,score]=pca(data.traces(idx,:));
            else
                score=data.score(idx,:);
                feats=params.feats;
            end
            clust_idx=kmeans(score(:,feats),clust_num); 
            clust_idx(clust_idx>1)=clust_idx(clust_idx>1)+max(clusters)-1;
            clust_idx(clust_idx==1)=split_clust; %leave one sub-clust with original number
            clusters(idx)=clust_idx;
        
            draw_clusters;
        end
    end

end
function outlier_split(hObj,event)
    global clusters handles data params;    

    axes(handles.axes1);
    colormatrix=colormap(lines);

    opt.WindowStyle='normal';
    answer=inputdlg({'Which Cluster to split outliers?','How many features?'},'Outlier removal',1,{'1','3'},opt);
    if(numel(answer) & numel(str2num(answer{1})))
        split_clust=str2num(answer{1});
        feat_num=str2num(answer{2});0
        idx=find(clusters==split_clust);
        if(numel(idx))
            mdist=sqrt(mahal(data.score(idx,1:feat_num),data.score(idx,1:feat_num))); %mahalanobis distance of samples
            h=figure;            
            [N,X]=hist(mdist,100);
            stem(X,N);            
            xlabel('Mahal. Distance');            
            ylabel('Frequency');
            axis([X(1) X(end) 0 10]);
            answer=inputdlg({'Outlier Threshold?'},'',1,{''},opt);
            if(numel(answer) & numel(str2num(answer{1})))
                ind=find(mdist>str2num(answer{1}));
                if(numel(ind))
                    clusters(idx(ind))=max(clusters)+1;
                    draw_clusters;
%                     draw_xcor;                
                end
            end
            close(h);
        end
    end

end

function merge_clusters(hObj,event)
    global  clusters handles data params;    
    axes(handles.axes1);
    colormatrix=colormap(lines);
    for i=1:max(clusters)
        names{i}=num2str(i);
    end
    [Selection,ok]=listdlg('ListString',names);
    if(ok)        
        for i=2:length(Selection);
            ind=find(clusters==Selection(i));
            clusters(ind)=Selection(1);
        end
        
        %compress clusters
        k=1;
        clust_idx=unique(clusters)
        for i=1:numel(clust_idx)
            c=clust_idx(i);
            if(c==0)
                continue;
            end
            ind=find(clusters==c);
            clusters(ind)=k;
            if(numel(ind))
                k=k+1;
            end            
        end

        draw_clusters;
%         draw_xcor;
    end
end
        
function show_shapes(hObj,event)
    global clusters handles data params;    
    grey=0.5*[1 1 1];

    Hf=figure;
    Ha=axes;
    colormatrix=colormap(lines);    
    maxtrace=0;
    mintrace=0;

    t=data.t;
            
    for i=1:max(clusters)               
        if(~numel(find(clusters==i)))
            continue;
        end

        mtrace=mean(data.traces(clusters==i,:),1);
        strace=std(data.traces(clusters==i,:),[],1);
        maxtrace=max(maxtrace,max(mtrace+strace));
        mintrace=min(mintrace,min(mtrace-strace));
        my_plotWithConf(t,mtrace,strace,colormatrix(i,:));
        hold on;
    end
    set(Ha,'XLim',[0 t(end)]);
    set(Ha,'YLim',[mintrace*1.05 maxtrace*1.3]);    
end

function remove_artifact(hObj,event)
%     clear clusters handles data params files;     
    global clusters handles data params;    

    for i=1:max(clusters)
       names{i}=num2str(i);
    end
    [Selection,ok]=listdlg('ListString',names); 
    if(ok)        
%             ind=1:size(data.traces,1);
            %move artifact times and traces to art data
        for i=1:length(Selection);
           s_ind=find(clusters==Selection(i));
           clusters(s_ind)=0;
        end
        
        %compress clusters
        k=1;
        clust_idx=unique(clusters)
        for i=1:numel(clust_idx)
            c=clust_idx(i);
            if(c==0)
                continue;
            end
            ind=find(clusters==c);
            clusters(ind)=k;
            if(numel(ind))
                k=k+1;
            end            
        end

        %re-compute PCA and draw
        display('Computing PCA.....');
        [data.coeff,data.score,latent]=pca(data.traces);
        data.eigvals=cumsum(latent)./sum(latent);
%         data.t=(1:size(data.traces,2))/params.samplerate;
%         axes(handles.axes1);
%         cla;
%         colormatrix=colormap(lines);    
        
%         plot(data.score(:,1),data.score(:,2),'.k');

        draw_clusters;
%         draw_xcor;
    end            
end

function export_data(hObj,event)
    global  clusters data outstr flag;        
    
%     clnums=unique(clusters);
    for i=1:max(clusters)
        clnames{i}=['Cluster ',num2str(i)];
        rastnames{i}=num2str(i);
    end
            
%     names=['Artifact',clnames];
    dlg_title='Event number for each cluster';
    numlines=1;
    defaultanswer=rastnames;
    answer=inputdlg(clnames,dlg_title,numlines,defaultanswer);            
        
    event_names=clusters;   
    for c=1:numel(answer);
        event_names(clusters==c)=str2num(answer{c});
        trace_m(c,:)=mean(data.traces(clusters==c,:),1);
        trace_s(c,:)=std(data.traces(clusters==c,:),[],1);
        isihist(c,:)=nanmean(data.isihist(clusters==c,:),1);
    end
    
    outstr.indices=event_names;
    outstr.traces.mean=trace_m;
    outstr.traces.std=trace_s;
    outstr.isihist=isihist;
    
    close(gcf);
    flag=1;
end

function savefile(file,strct,params)
    data=strct;
    save(file,'data','params');
end
%% utility functions

function draw_clusters()
    global clusters handles data params;    
    grey=0.5*[1 1 1];
    
    axes(handles.axes2);    
    colormatrix=colormap(lines);    
    cla;    
    t=data.t;
    maxtrace=max(data.traces(:));
    mintrace=min(data.traces(:));
            
    for i=1:max(clusters)               
        if(~numel(find(clusters==i)))
            continue;
        end
        
        H=plot(t,data.traces(clusters==i,:)');
        set(H,'Color',colormatrix(i,:));
        text('Position',[t(round(end/4)),maxtrace*1.1],'String',['cluster ',num2str(i)]...
            ,'Color',colormatrix(i,:),'FontSize',12);
        hold on;   
        t=t+size(data.traces,2)/params.samplerate;           
    end
    set(handles.axes2,'XLim',[0 t(1)]);
    set(handles.axes2,'YLim',[mintrace*1.05 maxtrace*1.3]);

    axes(handles.axes3);    
    cla;    
    bins=data.isibins;
    maxtrace=0;
    mintrace=0;
            
    for i=1:max(clusters)               
        if(~numel(find(clusters==i)))
            continue;
        end
        h=nanmean(data.isihist(clusters==i,:),1);
        H=plot(bins,h);
        maxtrace=max([maxtrace,max(h)]);
        set(H,'Color',colormatrix(i,:));
        hold on;   
        bins=bins+data.isibins(end);           
    end
    set(handles.axes3,'XLim',[0 bins(1)]);
    set(handles.axes3,'YLim',[mintrace*1.05 maxtrace*1.3]);
    
    axes(handles.axes1);    
    cla;
    colormatrix=colormap(lines);        
       
    for i=1:max(clusters)  
        iind=find(clusters==i);
        H=plot3(data.score(iind,params.feats(1)),data.score(iind,params.feats(2)),data.score(iind,params.feats(2)),'.');
        set(H,'Color',colormatrix(i,:));
        hold on;   
    end               
    
end

