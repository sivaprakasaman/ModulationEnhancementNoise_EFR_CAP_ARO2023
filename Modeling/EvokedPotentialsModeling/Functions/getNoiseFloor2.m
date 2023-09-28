function [floorx,floory] = getNoiseFloor2(pos_all,neg_all,Fs)
%Simple computation of noisefloor

num_trials = length(pos_all);
pos_trials = zeros(length(pos_all{1}),num_trials);
neg_trials = zeros(length(neg_all{1}),num_trials);
floorx = Fs*(0:length(pos_all{1})/2)/length(pos_all{1});
iters = 100;

%convert pos and neg cells to array

for i = 1:num_trials
    
    pos_trials(:,i) = pos_all{i};
    neg_trials(:,i) = neg_all{i};
    
end

floory_all = zeros(round(length(pos_all{1})/2),iters);

for j = 1:iters
    %randomly choose a subset of 1/5th number of trials:
    r_neg = randsample(num_trials,round(num_trials/6)); 
    r_pos = randsample(num_trials,round(num_trials/6)); 

    r_pos_trials = pos_trials(:,r_pos);
    r_neg_trials = neg_trials(:,r_neg);

    r_pos_trials(:,2:2:end) = -r_pos_trials(:,2:2:end);
    r_neg_trials(:,2:2:end) = -r_neg_trials(:,2:2:end);

    combined_avg = mean(horzcat(r_pos_trials, r_neg_trials),2);
    floor_fft = fft(combined_avg*1e6);

    floory_db = db(abs(floor_fft/length(pos_all{1})));
    floory_all(:,j) = floory_db(1:ceil(length(pos_all{1})/2));

end 

floory = mean(floory_all,2)';

end