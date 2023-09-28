function xx = transduc_nL(x,slope,asym)
    
    %x is vector here...better way to check this?
    %doesn't appear to depend on CF...even tho that's a param in the model
    
    corner = 80;
    strength = 20e6/(10^(corner/20));
    
    for i = 1:length(x)
        xx(i) = log(1.0+strength*abs(x(i)))*slope;
        if(x(i)<0)
            splx   = 20*log10(-x(i)/20e-6);
            asym_t = asym -(asym-1)/(1+exp(splx/5.0));
            xx(i) = -1/asym_t*xx(i);
        end
    end
    
end

