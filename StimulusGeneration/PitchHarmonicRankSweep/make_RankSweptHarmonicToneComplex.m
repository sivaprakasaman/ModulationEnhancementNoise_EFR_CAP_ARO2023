function [stim,time_sec] = make_RankSweptHarmonicToneComplex(amp,F0,fs,dur_sec,rate,N_harms,rel_amp,phi,start_rank,ramp)
    
    if ~exist('amp', 'var') || isempty(amp)
        amp = .95;
    end
    
    if ~exist('F0', 'var') || isempty(F0)
        F0 = 103;
    end
    
    if ~exist('fs', 'var') || isempty(fs)
        fs= 48828;
    end
    
    if ~exist('dur_sec', 'var') || isempty(dur_sec)
    	dur_sec = 1.3;
    end
    
    if ~exist('rate', 'var') || isempty(rate)
        rate = 12;
    end
    
    if ~exist('start_rank', 'var') || isempty(start_rank)
        start_rank = 1;
    end
    
    if ~exist('N_harms', 'var') || isempty(N_harms)
        N_harms = 6;
    end
    
    if ~exist('rel_amp','var') || isempty(rel_amp)
        rel_amp = ones(1,N_harms);
    end
    
    if ~exist('phi', 'var') || isempty(phi)
        phi = zeros(1,N_harms);
        phi(1:2:end) = pi/2; %keeping this the same as initial stim...maybe do 2:2:end??
    end
    
    if ~exist('ramp', 'var') || isempty(ramp)
       ramp = 0.02;
    end
    
    rate_hz = F0*rate;
    harm_no = (0:(N_harms-1))+start_rank;
    samples = 0:fs*dur_sec;
    
    time_sec = samples/fs;
    stim = zeros(1,length(samples));
    mags = rel_amp.*ones(1,N_harms);

%     mags(1) = mags(1)/2;
%     mags(end) = mags(end)/2;
    
    for i = 1:length(harm_no)
        
        stim = stim + mags(i).*sin(2*pi*F0*harm_no(i)*time_sec + pi*rate_hz*time_sec.^2 + phi(i));
        
    end

    stim = amp*stim/max(stim);
    stim = rampsound(stim,fs,ramp);
end