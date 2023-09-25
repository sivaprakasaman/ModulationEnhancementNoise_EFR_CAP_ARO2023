clear;
fc=503;
fm=103;
modDepth=1;
% duty=25;
dur=1;
Es=80;
fs=48828;
amp=20e-6*10.^(Es/20);

signal=make_SAM(fc, fm, fs, modDepth, dur, amp) ;
% signalrms=20e-6*10.^(Es/20);
% signal=(signalrms/(rms(signal))).*signal;



noise40=ltass_noise(fs,40,1,fs);
noise60=ltass_noise(fs,60,1,fs);

Stimulus_1=signal';
Stimulus_2=signal'+noise40;
Stimulus_3=signal'+noise60;
%% you shouldn't normalize !!! becusse signal level we want to be the same
% Stimulus_1=Stimulus_1/(1.1*max(Stimulus_1));
% Stimulus_2=Stimulus_2/(1.1*max(Stimulus_2));
% Stimulus_3=Stimulus_3/(1.1*max(Stimulus_3));

% %%
figure
subplot(1,3,1)
plot(Stimulus_1)
subplot(1,3,2)
plot(Stimulus_2)
subplot(1,3,3)
plot(Stimulus_3)
% %%
% rms(Stimulus_1)
% rms(Stimulus_2)
% rms(Stimulus_3)
audiowrite('SAM_503_0dB.wav',Stimulus_1,fs)
audiowrite('SAM_503_40dB.wav',Stimulus_2,fs)
audiowrite('SAM_503_60dB.wav',Stimulus_3,fs)

