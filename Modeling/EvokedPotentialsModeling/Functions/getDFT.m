function [f,DFT,PLV, floory] = getDFT(tot,collected,window,Fs,Fs0,gain,K_MRS, NF_iters)
%getDFT - Returns f - frequency, and mean DFT with noise floor removed.
%   Assumtions: 
%   -Pull 1/5th of total trials to look at and average 

len = length(tot)/2; 
numtrials = round(collected/5); %Number of trials to pull from/polarity, 1/5 of total collected

%% Separate out the +/- polarities

ind = 1;
pos = cell(1,len);
neg = cell(1,len);

%Change numtrials to l_SAM/sq25/etc if pooling randomly!
for i = 1:1:len
    %Pos 
    temp = tot{ind}(round(window(1)*Fs0):round(window(2)*Fs0))/gain;
    pos{i} = resample(temp,Fs,round(Fs0));
    
    %Neg
    temp2 = tot{ind+1}(round(window(1)*Fs0):round(window(2)*Fs0))/gain;
    neg{i} = resample(temp2,Fs,round(Fs0));
    
    ind = ind+2;
end

%fprintf('+/- Polarities Separated\n')

%% Calculate Mean Raw Spectrum


for i = 1:K_MRS
    
    [f,MRS(i,:),PLV(i,:)] = getSpectMag(pos,neg,Fs,numtrials); 
 %   fprintf('(Spectrum) Iteration %d of %d complete.\n',i,K_MRS)
    
end

MeanDFT = mean(MRS);
PLV = mean(PLV);

% plot(f,MeanDFT);

%% Calculate Noise Floor

[~, floory] = getNoiseFloor(pos,neg,Fs,NF_iters);

DFT = MeanDFT-floory;

end



