function generate_amp_file(amp,changroups,eodind,b_dhpf,a_dhpf)
%   GENERATE_AMP_FILE generate long *.amp file for FAST spike sorting
%   amp- input voltage matrix (column/channel)
%   changroups- cell array with vectors of channel ids for each tetrode to
%               be analyzed.
%   eodtimes- times of eod pulses

%% PARAMETERS
datapath=getenv('DATAPATH');
sdate=getenv('SESSDATE');
sesspath=[datapath,'\',sdate,'\'];
fname1=[sesspath,'\FAST\FAST_data'];
fid1 = fopen([fname1 '.amp'], 'a', 'l');
fname2=[sesspath,'FAST_idx'];

samplerate=str2num(getenv('SAMPLERATE'));
eodtimes=eodind/samplerate;
%blanking
blank_pre=2e-3*samplerate;
blank_post=8e-3*samplerate;

% high-pass 
% d = designfilt('highpassiir', ...       % Response type
%        'StopbandFrequency',100, ...     % Frequency constraints
%        'PassbandFrequency',150, ...
%        'StopbandAttenuation',55, ...    % Magnitude constraints
%        'PassbandRipple',4, ...
%        'DesignMethod','cheby1', ...     % Design method
%        'MatchExactly','stopband', ...   % Design method options
%        'SampleRate',samplerate)               % Sample rate
% [b_dhpf,a_dhpf]=tf(d);


%% 
K=10;
L=10;
I=eodind*ones(1,blank_pre+blank_post+1) + ones(size(eodind,1),1)*[-blank_pre:blank_post];
Io=eodind*ones(1,K*2) + ones(size(eodind,1),1)*[[(-blank_pre-K*L):L:(-blank_pre-1)] [(blank_post+1):L:(blank_post+K*L)]];
% Io=setdiff([1:size(amp,1)],I);
data=[];
Io=Io(:);
if(sum(Io>size(amp,1))>0)
    Io(Io>size(amp,1))=[];
    Io=[Io;size(amp,1)];
end
I=I(:);
if(sum(I>size(amp,1))>0)
    I(I>size(amp,1))=[];    
end

for i=1:numel(changroups)
    A=amp(:,changroups{i}); %take channels of current tetrode    
    for j=1:numel(changroups{i})
        a=A(:,j);
        %blank EODs
        B=interp1(Io(:),a(Io(:)),I,'pchip');
        a(I)=B;
        a=FiltFiltM(b_dhpf,a_dhpf,a);
        %blank EODs
%         a(I)=0;
        A(:,j)=a;
    end    

    data=[data A];
end

i=1;
nsamp=size(data,1);
nCh=size(data,2);
if(~exist(fname2))
    fid2 = fopen(fname2, 'a');
    fwrite(fid2,1,'uint64'); %first file- index=1
else
    fid2 = fopen(fname2, 'r');
    idx=fread(fid2,inf,'uint64'); %read all existing    
    fclose(fid2);
    fid2 = fopen(fname2, 'a');
    fwrite(fid2,idx(end)+nsamp,'uint64'); %index of this file within long file
end
fclose(fid2);

while(i<nsamp)
    clc;
    display(i/nsamp);
    I=min(nsamp,i+1e6);
    dat=data(i:I,:)'/1e6;
    % Pad data with zeros to 64 channels
    dat = cat(1, double(dat), zeros(64-nCh, size(dat,2)));

    % Convert to uint16 precision
    dat = uint16(dat/1.95e-7 + 32768);
    
    % Flatten data array
    chipInd = [1:32; 33:64];
    dat = dat(chipInd(:),:);
    dat = dat(:);

    fwrite(fid1, dat, 'uint16');
    i=I;
end
fclose(fid1);
end
