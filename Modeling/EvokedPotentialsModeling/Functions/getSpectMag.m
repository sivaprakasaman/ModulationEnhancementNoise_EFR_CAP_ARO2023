function [f,P1,PLV] = getSpectMag(pos,neg,Fs,numtrials)
%EDITED FOR MODELING!!!!!!!

%Returns the mean raw spectral magnitude in dB
% Example:
% pos = all positive trials
% neg = all negative trials
% Fs0 = 48828.125; %sampling rate in
% Fs = 15e3 (example); %resample to
% numtrials = 100; %Number of trials to pull per polarity

%% Parameters
len = size(pos,2); %number of trials collected/polarity

%% Pool random odds and evens from separated trials

%since already separated, just need to pick a random set of n trials from
%1-length of total #

%THIS WON"T WORK UNLESS ABOVE SECTION IMAX is changed!
r_odds = randi([1,len],[numtrials,1]);
r_evens = randi([1,len],[numtrials,1]);

% pos_r = zeros(length(pos{:,1}),numtrials);
% neg_r = zeros(length(neg{:,1}),numtrials);

for i = 1:1:numtrials

   pos_r(:,i) = pos(:,r_odds(i));
   neg_r(:,i) = neg(:,r_evens(i));
    
end

%% Remove DC

pos_r = pos_r-mean(pos_r);
neg_r = neg_r-mean(neg_r);

%% Create FFT Matrices of pos and neg data for analysis
%remember fft is applied to COLUMNS

pos_fft = fft(pos_r);
neg_fft = fft(neg_r);

sum_pos = sum(pos_fft,2);
sum_neg = sum(neg_fft,2);

sum_all = (sum_pos+sum_neg)/(2*numtrials);

L = length(sum_all);

P2 = 20*log10((abs(sum_all/(L)))); %Taking the average numtrials*2??
P1 = P2(1:floor(L/2)+1);

f = Fs*(0:(L/2))/L;


%% Calculate PLV
pos_phase = angle(pos_fft(1:floor(L/2)+1,:));
pos_vector = sum(exp(1i*pos_phase),2);

neg_phase = angle(neg_fft(1:floor(L/2)+1,:));
neg_vector = sum(exp(1i*neg_phase),2);

PLV = abs(pos_vector+neg_vector)/(numtrials*2);%maybe times2?

%% Calculate the mean of all trials
% %Setting this up in a way that should work well with the above logic
% 
% pos_sum = zeros([1,length(pos{1})]);
% neg_sum = zeros([1,length(neg{1})]);
% 
% %Compute the averages of pos and neg polarities, remove DC, and get a sum
% for i = 1:numtrials
%     pos_sum = pos_r{i} + pos_sum;
% end
% 
% mean_pos = (pos_sum-mean(pos_sum))/numtrials;
% 
% for i = 1:numtrials
%     neg_sum = neg_r{i} + neg_sum;
% end
% 
% mean_neg = (neg_sum-mean(neg_sum))/numtrials;
% 
% sum_all = mean_pos+mean_neg;
% %% FFT
% mag_envresponse = fft(sum_all*1e6);
% L = length(sum_all);
% P2 = 20*log10((abs(mag_envresponse/L)));
% P1 = P2(1:L/2+1);
% 
% f = Fs*(0:(L/2))/L;
% 
% %% PLV
% % 
%   pos_fft = fft(mean_pos*1e6);
%   L = length(sum_all);
%   pos_phase = angle(pos_fft(1:L/2+1));
%   pos_vector = exp(1i*pos_phase);
%   
%   neg_fft = fft(mean_neg*1e6);
%   L = length(sum_all);
%   neg_phase = angle(neg_fft(1:L/2+1));
%   neg_vector = exp(1i*neg_phase);
% 
% %PLV = angle(mag_envresponse);
% %PLV = PLV(1:L/2+1);

end

