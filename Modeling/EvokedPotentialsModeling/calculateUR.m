function [u_response,UR_fs,dur] = calculateUR(response_onset, ihc_loss_flag,CF,modelParams)

if ~ihc_loss_flag %normal
     load("a0012_ABR_click.mat");   % MP CHANGE THIS BASED ON WHAT AS GIVES YOU
     ABR = x.AD_Data.AD_Avg_V{1}./x.AD_Data.Gain;
     ABR_fs = round(x.Stimuli.RPsamprate_Hz);
else
    
end

fs = ABR_fs; % = 48828
t = (1:length(ABR))/ABR_fs;
%click_onset = 7; %7ms ABR response onset
amp_90dB = 20e-6 * 10^(90.0/20.0)/.7071;
input = zeros(size(t));
input(1,round(response_onset/1000*fs)) = 2*amp_90dB;
dur=max(t);
modelParams.dur = max(t);
%% Setting up Parameters

% n_cf=50;
% CF=logspace(log10(125.0), log10(20000.0), n_cf);

%Normal Model Params


%% Get APpsth for click

%feed x into model
[psth_pos,psth_neg,psth_fs] = getAP_PSTH_cihc(input,fs,modelParams,CF, ihc_loss_flag);
sum_psthpos = sum(psth_pos,1);
sum_psthneg = sum(psth_neg,1);
psth_sum = sum_psthpos+sum_psthneg;

%% Deconvolve

ABR_resamp = resample(ABR,psth_fs,ABR_fs);
% ABR_resamp = ABR_resamp./max(ABR_resamp);

psth_sum_sampled = psth_sum(1,1:length(ABR_resamp)); %match length of upsampled ABR
% psth_sum_sampled  = psth_sum_sampled./max(psth_sum_sampled);

u_response = deconvtvl2(ABR_resamp, psth_sum_sampled, 100);
UR_fs = psth_fs; 

end