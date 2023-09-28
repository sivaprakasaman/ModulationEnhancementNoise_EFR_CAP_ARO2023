function [sig, t] = make_RAM(fc, fm, fs, modDepth, duty, dur, amp, phi_c, phi_m) 
    if nargin <3
        error ('Need more love');
    end
    
    if ~exist('phi_m', 'var')
        phi_m= -pi/2;
    elseif isempty(phi_m)
        phi_m= -pi/2;
    end
    
    if ~exist('phi_c', 'var')
        phi_c= 0;
    elseif isempty(phi_c)
        phi_c= 0;
    end
    
    if ~exist('amp', 'var')
        amp= 1;
    elseif isempty(amp)
        amp= 1;
    end
    
    if ~exist('duty', 'var')
        duty= 25;
    elseif isempty(duty)
        duty= 25;
    end
    
    if ~exist('dur', 'var')
        dur= .2;
    elseif isempty(dur)
        dur= .2;
    end
    
    if ~exist('modDepth', 'var')
        modDepth= 1;
    elseif isempty(modDepth)
        modDepth= 1;
    end

    t = 0:1/fs:dur-1/fs;
    sig = amp*(modDepth*square(2*pi*fm*t + phi_m, duty)+2-modDepth)/2.*sin(fc*2*pi*t + phi_c);

end
