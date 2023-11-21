function [args_out] = make_EvokedPotential(stim,CFs,wts,UR,UR_fs, model_num)


%model_num
% 0 - efferent model
% 1 - BEZ2018
% 2 - ??

if ~exist("model_name","var")
    model_name = 0;
end

% add 2018 and a switch to use either models

switch model_name

    case 0
        disp('Simulating  using SIM_EFFERENT_MODEL');
        [p_ihcout_s, p_hsr_s, p_lsr_s, p_ic_s, p_gain_s] = sim_efferent_model(stim,CFs);
        [n_ihcout_s, n_hsr_s, n_lsr_s, n_ic_s, n_gain_s] = sim_efferent_model(-stim,CFs);
    case 1
        disp('Simulating using BEZ2018 model. No Efferent, no IC.');
            
        %TODO:
            % - Ensure sim_efferent_model uses similar setup params or vice
            % versa
            % - 

        %KEY Differences:
            % - No IC response, only simulating an HSR fiber (can modify if
            % needed)
            % - Where to calculate level??

            dB_stim = 75;
            modelParams.tabs = 0.6e-3;
            modelParams.trel = 0.6e-3;
            modelParams.cohc = ones(length(CF),1); %healthy
            modelParams.cihc = ones(length(CF),1); %healthy
            modelParams.species = 1; % 1 for cat (2 for human with Shera et al. tuning; 3 for human with Glasberg & Moore tuning) %read up on this tuning
            modelParams.noiseType = 0; % 1 for variable fGn; 0 for fixed (frozen) fGn
            modelParams.implnt = 0; % "0" for approximate or "1" for actual implementation of the power-law functions in the Synapse
            modelParams.stimdb = dB_stim;
            modelParams.reps = 100; 
            modelParams.Fs = 100e3;
            modelParams.psthbinwidth = 1e-4;
            modelParams.buffer = 2;
            modelParams.cihc = 1;
            modelParams.cohc = 1;
            modelParams.spont = 55;
            



end


% we should explore this window
win=[0.1,0.9];

%reshape to a new sample rate, only considering HSRs for now.
ap_psth_env_full = (p_hsr_s+n_hsr_s)/2;
ap_psth_env_full=ap_psth_env_full';
ap_psth_env_ic_full = (p_ic_s+n_ic_s)/2;
ap_psth_env_ic_full=ap_psth_env_ic_full';

% we should define ap_psth_env and ap_psth_env_ic

psthbins=10;
for i=1:length(CFs)
    ap_psth_env(:,i) = sum(reshape(ap_psth_env_full(:,i),psthbins,length(ap_psth_env_full(:,i))/psthbins));
    ap_psth_env_ic(:,i)= sum(reshape(ap_psth_env_ic_full(:,i),psthbins,length(ap_psth_env_ic_full(:,i))/psthbins));
end

fs_efr = 100e3/psthbins; %model sample rate/bins

%normalized
an_cap = sum(ap_psth_env(win(1)*fs_efr:win(2)*fs_efr,:),2);
an_cap = an_cap/max(an_cap);
an_cap = an_cap-mean(an_cap);


ic_cap = sum(ap_psth_env_ic(win(1)*fs_efr:win(2)*fs_efr,:),2);
ic_cap = ic_cap/max(ic_cap);
ic_cap = ic_cap-mean(ic_cap);

%TODO: convolve with click later???

comb_cap = wts(1).*an_cap + wts(2).*ic_cap;
UR=resample(UR,fs_efr,UR_fs);
UR=UR-mean(UR);
convout=conv(comb_cap, UR);

args_out.convout = convout;
args_out.ic_cap = ic_cap;
args_out.an_cap = an_cap;
args_out.comb_cap = comb_cap;
args_out.fs_cap = fs_efr;

end