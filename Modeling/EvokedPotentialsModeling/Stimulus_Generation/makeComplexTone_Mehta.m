function x = makeComplexTone_Mehta(F0, dur, fs, db_main, db_flank, rank, nharms, ramp, phi, noise_on)

if(~exist('fs','var'))
    fs = 48828.125; % Sampling Rate
end

if(~exist('db_main','var'))
    db_main = 55; % dB of 10 middle tones
end

if(~exist('db_flank','var'))
    db_flank = 49; % dB of 2 flanked tones
end

if(~exist('rank','var'))
    rank = 10; % dB of 2 flanked tones
end

if(~exist('nharms','var'))
    nharms = 10;
end

if(~exist('dur','var'))
    dur = 1; % Duration in Seconds
end

if(~exist('ramp','var'))
    ramp = 0.010; %In seconds
end

if(~exist('F0','var'))
    F0 = 440;
end

if(~exist('phi','var'))
    phi = ones(nharms+2, 1) * pi/2; % Sin for all harmonics
end

if(~exist('noise_on','var'))
    noise_on = 0; % noise off by default
end

t = 0:(1/fs):(dur - 1/fs);
x = 0;

harmonics = (rank-1):(rank+nharms);

for k = 1:length(harmonics)
    
    if(k==1 || k==length(harmonics))
        mag = db2mag(db_flank); 
    else
        mag = db2mag(db_main);
    end
    
    x = x + mag*cos(2*pi*F0*harmonics(k)*t + phi(k));

end

x = scaleSound(rampsound(x, fs, ramp));

end

