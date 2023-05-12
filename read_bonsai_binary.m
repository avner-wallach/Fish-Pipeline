function [ t,data ] = read_bonsai_binary(fileName,samplefreq,numChannels,blocksize,chans,source)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
T=31*60; %time of recording (in s)
N=50; %number of blocks per channel in each read
if(all(chans<2)) %indices
    chans=find(chans);
end
if(nargin<6)
    chans=1:numChannels;
end
if(nargin<5)
    blocksize=256;
end
bits=16;

%initialize
A=uint16(zeros(T*samplefreq,numel(chans)));

%read it
D = fopen(fileName);
% tic;
B = fread(D,N*blocksize*numChannels,'uint16');
k=1; j=1;
while(numel(B)>0)    
    B=reshape(B,blocksize,[]);
    %truncate in case saving crashed in mid saving
    b=size(B,2);
    B=B(:,1:(b-mod(b,numChannels)));
    M=size(B,2)/numChannels*size(B,1);    
    for i=1:numel(chans)
        try
            A(k:(k+M-1),i)=reshape(B(:,chans(i):numChannels:end),[],1);
        catch
            display(size(B));
        end
    end
    k=k+M;    
%     clc;
%     disp(k/(T*samplefreq));   
%     a=toc/N;
%     disp(a)
%     tic
    B = fread(D,N*blocksize*numChannels,'uint16');
end

A(k:end,:)=[]; %remove trailing zeros

% j=1;
% for k=1:numel(C)
%     if(ismember(k,chans))
%         A(:,j)=C{k};
%         j=j+1;
%     end
% end
% B = fread(D,Inf,'uint16'); %read in m x n matrix, m = numChannels
% B=reshape(B,blocksize,[]);
% %truncate in case saving crashed in mid saving
% b=size(B,2);
% B=B(:,1:(b-mod(b,numChannels)));
% for i=1:numel(chans)
%     A(:,i)=reshape(B(:,chans(i):numChannels:end),[],1);
% end
    
% %sanity
% if size(A,2) ~= numChannels
%     error('numChanels does not match desired raw data import')
% else
%     warning('sanity check input passed')
% end
S=2^(bits-1);
if(strcmp(source,'aux'))
    a=0.0000374;%in v
elseif(strcmp(source,'adc'))
    a=0.195; %in uV
else
    a=0.1007;
end
data=(double(A)-S)*a;

t=[1:size(A,1)]/samplefreq;

fclose(D);
% 
% %plot it
% figure,
% kk_stackedLinePlot(A(1:numChannels/4,1:1e4)), title(sprintf('plot of %s',fileName)); %plot first 8 channels, first 1e4 samps
% xlabel('samples'), ylabel('Chan#')
end

