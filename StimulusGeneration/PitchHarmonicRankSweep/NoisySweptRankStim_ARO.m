
%how much quieter in dB the noise should be relative to signal.
db_drop = 10; 

fs = 48828;
amp = 0.1;
rms_sig = amp;
F0 = 103;
dur_sec = 1;
rate = 12;
N_harms = 6;
rel_amp = ones(1,N_harms);

%ALT Phase
phi = zeros(1,N_harms);
phi(1:2:end) = pi/2; 

start_rank = 1;
ramp = 0.02;

[sig, time_sec] = make_RankSweptHarmonicToneComplex(amp,F0,fs,dur_sec,rate,N_harms,rel_amp,phi,start_rank,ramp);
noise = ltass_noise(fs, 40, 1, dur_sec*fs+1);

sig = rms_sig*sig/rms(sig);
rms_noise = db2mag(mag2db(rms_sig)-db_drop);
noise = rms_noise*noise/rms(noise);

if isnan(db_drop)
    noise = 0;
    db_drop = [];
end

sig_out = sig+noise;
rev_sig_out = flip(sig_out);

audiowrite(['SweptRank_F0_', num2str(F0),'_ltassNoise_dbdrop_',num2str(db_drop),'_rate_',num2str(rate),'.wav'],sig_out,fs);
audiowrite(['REV_SweptRank_F0_', num2str(F0),'_ltassNoise_dbdrop_',num2str(db_drop),'_rate_',num2str(rate),'.wav'],rev_sig_out,fs);
