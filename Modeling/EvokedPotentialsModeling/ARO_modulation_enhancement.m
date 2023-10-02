clear
close all

signallevels=40:10:90;
SNRlevels=-10:10:30;


for i=1:length(signallevels)
    for j=1:length(SNRlevels)
      
        signallevel=signallevels(i);
        SNR=SNRlevels(j);
        
        fc=503;
        fm=103;
        modDepth=1;
        % duty=25;
        dur=1;
        Es=signallevel;
        fs=1e5;
        amp=20e-6*10.^(Es/20);

        signal=make_SAM(fc, fm, fs, modDepth, dur, amp) ;
       
        noise=ltass_noise(fs,signallevel-SNR,1,fs);
      

        % fix make_SAM so no neend for inversion (RAM doesn;t need inversion)
        SAM_sig=signal';
        SAM_signoise=signal'+noise;
      

        % we should discusse if this is the best range to use
        n_cf=50;
        CFs=logspace(log10(101.0), log10(20000.0), n_cf);

       
        %AN and IC weights, make sure sum to 1, this should be explored
        wts = [.5,.5];

        % We should try human ABR and maybe CA ABR and also level
        load("a0012_ABR_click.mat");   % normal chin ABR 80dB?
        ABR = x.AD_Data.AD_Avg_V{1};
        ABR_fs = round(x.Stimuli.RPsamprate_Hz);

        args_out_sig = make_EvokedPotential(SAM_sig,CFs,wts,ABR,ABR_fs);
        args_out_signoise = make_EvokedPotential(SAM_signoise,CFs,wts,ABR,ABR_fs);

        % to do add convout it creats error becuase size is different
        S=  [args_out_sig.comb_cap,args_out_sig.ic_cap,args_out_sig.an_cap];
        N = [args_out_signoise.comb_cap,args_out_signoise.ic_cap,args_out_signoise.an_cap];

        %Compute FFT for each comb_cap, ic_cap, and an_cap, sig and noise

        %find whether the noise is enhanced in the spectrum for each one


        nfft = 2^nextpow2(size(S,1));
        L = nfft;
        f = linspace(0,args_out_sig.fs_cap/2,nfft/2);

        Y_S = fft(S,nfft);
        Y_N = fft(N,nfft);

        %         f = 8000*(0:(L/2))/L;
        P2_S = abs(Y_S/L);
        P1_S = P2_S(1:L/2,:);
        P1_S(2:end-1,:) = 2*P1_S(2:end-1,:);

        P2_N = abs(Y_N/L);
        P1_N = P2_N(1:L/2,:);
        P1_N(2:end-1,:) = 2*P1_N(2:end-1,:);

        % change this to length of S and N
        for u = 1:3
            [pks_S(:,u),~,~] = getPeaks(f,P1_S(:,u),fm,1);
            [pks_N(:,u),~,~] = getPeaks(f,P1_N(:,u),fm,1);
        end
         % we should define this in the begining of the code
        modulationEnhancement(i,j,:) = pks_N>pks_S;

    end
end
   save('modulationenhancement.mat','modulationEnhancement')
