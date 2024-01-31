clear
close all

cwd = pwd();
%where to save simulations
simdir = 'SimulatedData';

if ~isfolder(simdir)
    mkdir(simdir);
end


% we should discusse if this is the best range to use
n_cf=25;
CFs=logspace(log10(125.0), log10(20000.0), n_cf);
CF=CFs;

modelParams.tabs = 0.6e-3;
modelParams.trel = 0.6e-3;
modelParams.cohc = ones(length(CFs),1); %healthy
modelParams.cihc = ones(length(CFs),1); %healthy
modelParams.species = 1; % 1 for cat (2 for human with Shera et al. tuning; 3 for human with Glasberg & Moore tuning) %read up on this tuning
modelParams.noiseType = 0; % 1 for variable fGn; 0 for fixed (frozen) fGn
modelParams.implnt = 0; % "0" for approximate or "1" for actual implementation of the power-law functions in the Synapse
%modelParams.stimdb = dB_stim;
modelParams.reps = 50;
modelParams.Fs = 100e3;
modelParams.psthbinwidth = 1e-4;
modelParams.buffer = 2;
modelParams.spont = 55; %HSR

ihc_loss_flag = 0;
[unitary_resp, UR_fs,~] = calculateUR(7,ihc_loss_flag,CFs,modelParams);
plot(unitary_resp.f)

stim_dur = 1;
modelParams.dur = stim_dur;  %% ??

%%


ihc_grades = [100, 50,10,3]/100;

load("pitch_harmonic_rank_alt.mat");
Fs = 100e3;
fs=fs;
fc = 4000;
fm = 103;
modDepth = 1;
%dur = length(pitches(:,1))/fs; %use same dur
signallevel = 90;
Es=signallevel;
amp=20e-6*10.^(Es/20);
stim_dur = 1;

[sam_tone, t]= make_SAM(fc, fm, fs, modDepth, stim_dur, amp,[],[]);
[r25, t] = make_RAM(fc, fm, fs, modDepth, 25, stim_dur, amp, [], []);
[r50, t] = make_RAM(fc, fm, fs, modDepth, 50, stim_dur, amp, [], []);
%
%should have all my stims
all_stims = [pitches,sam_tone,r50',r25'];

% args_out_sig = make_EvokedPotential_new(sam_tone,CFs,unitary_resp,UR_fs,modelParams,ihc_loss_flag); % add ABR_fs as o/p, CFs as i/p

%%
for r = 1:size(all_stims,2)
    
    modelParams.cihc = ones(length(CF),1); %healthy
    fprintf('\n STIMULUS %i of %i',r,size(all_stims,2));
    input = all_stims(:,r);
    input = input';
    modelParams.dur = length(input)/fs;
    fprintf('\n    Normal Run');
    args_out_sig = make_EvokedPotential_new(input,CFs,unitary_resp,UR_fs,modelParams,ihc_loss_flag); % add ABR_fs as o/p, CFs as i/p
    grand_envs_n(:,r) = args_out_sig;
    
    for g = 1:length(ihc_grades)
        fprintf('\n    Impaired Run %i of %i', g, length(ihc_grades));
        modelParams.cihc = ihc_grades(g).*ones(length(CF),1);
        args_out_sig = make_EvokedPotential_new(input,CFs,unitary_resp,UR_fs,modelParams,1); % add ABR_fs as o/p, CFs as i
        grand_envs_i(g,:,r) =  args_out_sig;
    end
    
end

timestr = [num2str(yyyymmdd(datetime)),num2str(hour(datetime)),num2str(minute(datetime))];
cd(simdir);
save(['sim_',timestr,'.mat']);
cd(cwd);

