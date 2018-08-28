function [ vout ] = remove_60hz(t,vin)
% REMOVE_60HZ actively removes 60 hz noise by fitting 

%LPF
samplerate=1/mean(diff(t));
dnotch = designfilt('bandstopiir','FilterOrder',2, ...
    'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
    'DesignMethod','butter','SampleRate',samplerate);

vin1=filtfilt(dnotch,vin);

dlp = designfilt('lowpassiir','FilterOrder',8, ...
         'PassbandFrequency',2e3,'PassbandRipple',0.2, ...
         'SampleRate',samplerate);

 vin2=filtfilt(dlp,vin1);
end

