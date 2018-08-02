function [ t,data ] = read_bonsai_binary(fileName,samplefreq,numChannels,blocksize,chans,source)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
if(nargin<6)
    chans=1:numChannels;
end
if(nargin<5)
    blocksize=256;
end
bits=16;
%read it
% fileName = 'amplifier2014-11-25T22_09_28.bin'; %amplifierData, numChannels = 32;
D = fopen(fileName);
% numChannels = 32;
B = fread(D,Inf,'uint16'); %read in m x n matrix, m = numChannels
B=reshape(B,blocksize,[]);
for i=1:numel(chans)
    A(:,i)=reshape(B(:,chans(i):numChannels:end),[],1);
end
    
% %sanity
% if size(A,2) ~= numChannels
%     error('numChanels does not match desired raw data import')
% else
%     warning('sanity check input passed')
% end
S=2^(bits-1);
if(strcmp(source,'aux'))
    a=0.0000374;%in v
else
    a=0.195; %in uV
end
data=(A-S)*a;

t=[1:size(A,1)]/samplefreq;

fclose(D);
% 
% %plot it
% figure,
% kk_stackedLinePlot(A(1:numChannels/4,1:1e4)), title(sprintf('plot of %s',fileName)); %plot first 8 channels, first 1e4 samps
% xlabel('samples'), ylabel('Chan#')
end

