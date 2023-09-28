clear all 
close all

%Generate a set of stimuli:

F0 = 103;
ranks = [3,5,7,9,11,13,15,20,25];
nharms = 4;
phase = 'sin';
total_dur = 1;
dur = .30;
fs = 44e3;
db_main = 52;
db_flank = 46;
ramp = 0.02;

if strcmp(phase,'sin')
    phi = ones(nharms+2, 1) * pi/2;
elseif strcmp(phase,'alt')
    phi = zeros(nharms+2,1);
    phi(2:2:end) = pi/2;
end

figure;
hold on;
for i = 1:length(ranks)
    x = makeComplexTone_Mehta(F0, dur, fs, db_main, db_flank, ranks(i), nharms, ramp, phi);
    [Pxx, f] = pmtm(x,[],[],fs);
    plot(f,Pxx);
    pitches(:,i) = x;
end
hold off;

save(['pitch_harmonic_rank_',phase,'.mat'],'F0','fs','pitches','ranks')

