function []=save_raw_spfile(amp,ops,segnum,eodind,fname)
%   GET_SPIKE_TIMES extract spike times from continuous tetrode recording.
%   t- input time vector
%   amp- input voltage matrix (column/channel)
%   changroups- cell array with vectors of channel ids for each tetrode to
%               be analyzed.
%   eodtimes- times of eod pulses
%% PARAMETERS
samplerate=ops.samplerate;

%EOD blanking
eod_pre=ops.eodblankpre*samplerate;
eod_post=ops.eodblankpost*samplerate;

%LFP blanking
blank_pre=ops.blankpre*samplerate;
blank_post=ops.blankpost*samplerate;

%changroups
changroups=cell2mat(ops.seg(segnum).LFPgroups(ops.seg(segnum).Spkgroups));
A=amp(:,changroups); %take channels of current tetrode    

%filename
filename=[ops.spdatpath,'\spdata_seg',num2str(segnum),'.bin'];
logfile=[ops.spdatpath,'\spdata_seg',num2str(segnum),'.log'];
%% old method filtering (interpolation + blanking)
% K=100;
% 
% %blanking sections
% I=eodind*ones(1,blank_pre+blank_post+1) + ones(size(eodind,1),1)*[-blank_pre:blank_post];
% 
% %100 samples before and after sections
% Io=eodind*ones(1,K*2) + ones(size(eodind,1),1)*[[(-blank_pre-K):(-blank_pre-1)] [(blank_post+1):(blank_post+K)]];
% 
% I=unique(I(:));
% Io=unique(Io(:));
% Io=setdiff(Io,intersect(I,Io));
% I=I(I>0 & I<size(amp,1));
% Io=Io(Io>0 & Io<size(amp,1));
% if(I(1)==1)
%     Io=[1;Io];
% end
% for j=1:numel(changroups)
%     a=A(:,j);
%     %blank EODs
%     B=interp1(Io(:),a(Io(:)),I,'pchip');
%     a(I)=B;
%     A(:,j)=a;
% end    
%% new method filtering : median filtering + EOD blanking + LFP template cancelation
K=2e-3*samplerate; %med filt kernel
A=A-medfilt1(A,K,[],1); %use median filter to remove slow components

%EOD blanking sections
I_eod=eodind*ones(1,eod_pre+eod_post+1) + ones(size(eodind,1),1)*[-eod_pre:eod_post];
I_eod=unique(I_eod(:));
I_eod=I_eod(I_eod>0 & I_eod<size(amp,1));

A(I_eod(:),:)=0; %blank eod

%LFP blanking sections
eodind=eodind(eodind>blank_pre & eodind<(size(amp,1)-blank_post));
I=eodind*ones(1,blank_pre+blank_post+1) + ones(size(eodind,1),1)*[-blank_pre:blank_post];

for i=1:size(A,2)
    a=A(:,i);
    Y=a(I)';    
    [x,tr]=plot_graded(1:size(I,2),Y',[],max(abs(Y)),50,0);
    cc=corr(Y,tr');
    [c0,mm]=max(cc,[],2);
    T=tr(mm,:)';
    beta=(sum(T.*Y))./vecnorm(T).^2;
    D=beta'*ones(1,size(tr,2)).*T';    
%     medtrace(:,i)=nanmedian(a(I),1);
%     cr(i,:)=medtrace(:,i)\a(I)';
%     D=cr(i,:)'*medtrace(:,i)';
    a(I)=a(I)-D;
    A(:,i)=a;    
end

%% save data
% add empty channels if group is smaller than 4 channels
if(size(A,2)<4)
    A=[A zeros(size(A,1),4-size(A,2))];
end

%append data to binary file
F=fopen(filename,'a');
fwrite(F,A','int16');
fclose(F);

%append data to log file
% F=fopen(logfile,'a');
% fwrite(F,fname,'int16');
dlmwrite(logfile,fname,'-append');
% fclose(F);

end