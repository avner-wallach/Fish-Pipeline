function [ eodrast ] = get_eod_rasters( data,numunits )
%GET_EOD_RASTERS generate eod response raster from sorted data
rasterpre=str2num(getenv('RASTERTPRE'));
rasterpost=str2num(getenv('RASTERTPOST'));

% samplerate=str2num(getenv('SAMPLERATE'));
%% initialize
S=numel(data.SPIKES);
offset=data.FILE.offset;
for s=1:S
%     U=max(data.SPIKES(s).raster(:,2));
    U=numunits(s);
    eodrast{s}=cell(U,1);
    for u=1:U
        eodrast{s}{u}=zeros(sum(data.SPIKES(1).raster(:,2)==u)*10,2);
    end
end
% spcount=cell(S,1);
%%
teod=data.EOD.t;
T=numel(teod);
for s=1:S   %every electrode group
    raster=data.SPIKES(s).raster;    
    U=numunits(s);    
%     U=max(raster(:,2)); %number of units
%     spcount{s}=zeros(numel(teod),U);
    k=ones(1,U);
    for t=1:T  %every eod
        tt=teod(t);%+offset;
        rast=raster(inrange(raster(:,1),tt+[-rasterpre rasterpost]),:);
        for u=1:U
            ind=find(rast(:,2)==u);
            I=numel(ind);
            eodrast{s}{u}(k(u)+[0:(I-1)],:)=[rast(ind,1)-tt ones(numel(ind),1)*t];
            k(u)=k(u)+I;
%             spcount{s}(t,u)=sum(inrange(rast(ind,1)-tt,scwin));
        end
    end
    for u=1:U
        eodrast{s}{u}(k(u):end,:)=[];
    end
end

end

