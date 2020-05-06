function [fish,seg,coornames,data]=get_posture(txt,num,position)
% get DeepLabCut tracking data and convert to model
    medfiltk=str2num(getenv('MEDFILTK')); %tracking median kernel

    fnames=unique(txt(2,2:end));
    for i=1:numel(fnames)
        ind=find(strcmp(txt(2,:),fnames{i}));    
%         if(~strcmp(txt(1,1),'bodyparts'))
%             ind=ind+1;
%         end
        if(numel(ind))
            data.(fnames{i}).xy=num(:,(ind(1):ind(2)));
            data.(fnames{i}).c=num(:,(ind(3)));
        end
        if(nargin>2 & numel(position)) %video was cropped
            data.(fnames{i}).xy=data.(fnames{i}).xy+position;
        end
    end
    
    %% get midpoint between LED and Trunk1
    data.midpoint.xy=[mean([data.mouth.xy(:,1),data.Trunk1.xy(:,1)],2) ...
        mean([data.mouth.xy(:,2),data.Trunk1.xy(:,2)],2)];
    data.midpoint.c=min([data.mouth.c data.Trunk1.c],[],2);
    
    %% coordinate names
    coornames={'X','Y','azim','b_chin','b_tr2','b_tl1','b_tl2','b_cf',...
        'c_rt','c_lt'};

    th=0.8;
    %% extract location (I leave Z position for later)
    %X
%     fish(:,strcmp(coornames,'X'))=data.LED.xy(:,1);
%     fish(data.LED.c<th,strcmp(coornames,'X'))=NaN;
%     xy=(data.mouth.xy+data.Trunk1.xy)/2;
%     fish(:,strcmp(coornames,'X'))=xy(:,1);
%     fish(data.Trunk1.c<th | data.mouth.c<th,strcmp(coornames,'X'))=NaN;
    fish(:,strcmp(coornames,'X'))=data.midpoint.xy(:,1);
    fish(data.midpoint.c<th,strcmp(coornames,'X'))=NaN;

    %Y
%     fish(:,strcmp(coornames,'Y'))=data.LED.xy(:,2);
%     fish(data.LED.c<th,strcmp(coornames,'X'))=NaN;
%     fish(:,strcmp(coornames,'Y'))=xy(:,2);
%     fish(data.Trunk1.c<th | data.mouth.c<th,strcmp(coornames,'Y'))=NaN;
    fish(:,strcmp(coornames,'Y'))=data.midpoint.xy(:,2);
    fish(data.midpoint.c<th,strcmp(coornames,'Y'))=NaN;

%     center.xy=xy;
%     center.c=min([data.Trunk1.c,data.mouth.c]);
    
    %azim
%     mouthv=data.mouth.xy-data.LED.xy;
    mouthv=data.mouth.xy-data.Trunk1.xy;
    azim=atan2(mouthv(:,2),mouthv(:,1));
    fish(:,strcmp(coornames,'azim'))=azim;
    fish(data.Trunk1.c<th | data.mouth.c<th,strcmp(coornames,'azim'))=NaN;
%     fish(data.LED.c<th | data.mouth.c<th,strcmp(coornames,'azim'))=NaN;

    
    %beta chin
    fish(:,strcmp(coornames,'b_chin'))=get_angle(data.Trunk1,data.mouth,data.chin);

%     %beta trunk1
%     fish(:,strcmp(coornames,'b_tr1'))=get_angle(data.mouth,data.LED,data.Trunk1);

    %beta trunk2
    fish(:,strcmp(coornames,'b_tr2'))=get_angle(data.mouth,data.Trunk1,data.Trunk2);

    %beta tail1
    fish(:,strcmp(coornames,'b_tl1'))=get_angle(data.Trunk1,data.Trunk2,data.Tail1);

    %beta tail2
    fish(:,strcmp(coornames,'b_tl2'))=get_angle(data.Trunk2,data.Tail1,data.Tail2);

    %beta cuadal fork
    fish(:,strcmp(coornames,'b_cf'))=get_angle(data.Tail1,data.Tail2,data.CaudalFork);

%     %gamma right pect base
%     fish(:,strcmp(coornames,'c_rb'))=get_angle(data.LED,data.midpoint,data.RPecBase);

    %gamma right pect tip
    fish(:,strcmp(coornames,'c_rt'))=get_angle(data.mouth,data.midpoint,data.RPecTip,'r');

%     %gamma left pect base
%     fish(:,strcmp(coornames,'c_lb'))=get_angle(data.LED,data.midpoint,data.LPecBase);

    %gamma left pect tip
    fish(:,strcmp(coornames,'c_lt'))=get_angle(data.mouth,data.midpoint,data.LPecTip,'l');
    
    %% extract segment length
    seg(1)=get_segment(data.Trunk1,data.mouth)/2; %LED-mouth
    seg(2)=get_segment(data.mouth,data.chin); %mouth-chin
    seg(3)=get_segment(data.mouth,data.Trunk1)/2; %LED-trunk1
    seg(4)=get_segment(data.Trunk1,data.Trunk2); %trunk1-trunk2
    seg(5)=get_segment(data.Trunk2,data.Tail1); %trunk2-tail1
    seg(6)=get_segment(data.Tail1,data.Tail2); %tail1-tail2
    seg(7)=get_segment(data.Tail2,data.CaudalFork); %tail2-cf
    seg(8)=get_segment(data.midpoint,data.RPecTip); %center-rpt
    seg(9)=get_segment(data.midpoint,data.LPecTip); %center-lpt    
%     seg(8)=get_segment(data.LED,data.RPecBase); %LED-rb
%     seg(9)=get_segment(data.RPecBase,data.RPecTip); %rpb-rpt
%     seg(10)=get_segment(data.LED,data.LPecBase); %LED-lb
%     seg(11)=get_segment(data.LPecBase,data.LPecTip); %lpb-lpt

    function theta=get_angle(fielda,fieldb,fieldc,opt)
        if(nargin<4)
            opt=' ';
        end
        A=fielda.xy;
        B=fieldb.xy;
        C=fieldc.xy;    
        ca=fielda.c;
        cb=fieldb.c;
        cc=fieldc.c;
        AB=B-A;
        BC=C-B;
        theta=atan2(BC(:,2),BC(:,1))-atan2(AB(:,2),AB(:,1));
        if(opt=='l')
            theta=pi/2-theta;
        elseif(opt=='r')
            theta=theta+pi/2;
        end
        theta=mod(theta,2*pi);
        theta=mod(theta+pi,2*pi)-pi;        
        theta(ca<th | cb<th | cc<th)=NaN;
        theta(abs(theta)>pi/2)=NaN; %remove sporious results
        theta=medfilt1(theta,medfiltk,[],1,'omitnan');
    end

    function N=norm2(vec)
        N=hypot(vec(:,1),vec(:,2));
    end
    
    function seg=get_segment(fielda,fieldb)
        A=fielda.xy;
        B=fieldb.xy;
        ca=fielda.c;
        cb=fieldb.c;
        AB=B-A;
        N=norm2(AB);
        N(ca<th | cb<th)=NaN;
        seg=nanmean(N);
    end
        
    
end