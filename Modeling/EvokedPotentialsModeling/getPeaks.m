function [PKS,LOCS,n_floor] = getPeaks(f,DFT,F0,harmonics)
%getSum Outputs the sum of first 5 harmonics of a given .
%deviation. Also makes a plot.
%f = frequency vector
%DFT = DFT.
%F0 = fundamental freq
%harmonics = number of harmonics to pull 

% [PKS, LOCS] = findpeaks(DFT,f,'MinPeakHeight',max(DFT)/3,'MinPeakDistance',80);
% PKS = PKS.*(LOCS>99);
% LOCS = LOCS(1:harmonics);
% returns flat noise floor estimate based on non peak values 
% (use with caution!!!!!, crappy way of computing the floor)

%%harms = F0:F0:harmonics*F0;
PKS = zeros(1,harmonics);
LOCS = zeros(1,harmonics);
for_floor_inds = 1:length(DFT);

for i=1:harmonics
    
    %figure out the +/- range of frequencies to take max
    %to do:  10 should be a parameter input 
    r = DFT'.*((f>(i*F0-10)).*(f<(i*F0+10)));
    [PKS(i),ind] = max(r);
    LOCS(i) = ind;
%     exclude = ind-2:ind+2;
%     for j = 1:length(exclude)
%         for_floor_inds = for_floor_inds(for_floor_inds~=exclude(j));
%     end
    for_floor_inds = for_floor_inds(for_floor_inds~=ind);
    
end

%added for modeling to exclude crazy low freq noise floor variability
for_floor_inds = for_floor_inds(f(for_floor_inds)>F0);

n_floor = mean(DFT(for_floor_inds)).*ones(1,length(f));
% n_floor = spline(f(for_floor_inds),DFT(for_floor_inds),f);
harmsum = cumsum(PKS(1:harmonics));
PKS = PKS(1:harmonics);
harmsum = harmsum(1:harmonics);

%plot
%findpeaks(DFT,f,'MinPeakHeight',11,'MinPeakDistance',80)

end