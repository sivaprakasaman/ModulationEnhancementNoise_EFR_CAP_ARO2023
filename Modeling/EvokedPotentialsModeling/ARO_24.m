clear
close all

cwd = pwd();
%where to save simulations
simdir = 'SimulatedData';

if ~isfolder(simdir)
    mkdir(simdir);
end


% we should discusse if this is the best range to use
n_cf=50;
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
[unitary_resp, ABR_fs,dur] = calculateUR(7,0,CFs,modelParams);
plot(unitary_resp.f)
% args_out_sig = make_EvokedPotential_new(SAM_sig,CFs,unitary_resp,ABR_fs,modelParams,ihc_loss_flag); % add ABR_fs as o/p, CFs as i/p
% to do add convout it creats error becuase size is different
%S=  [args_out_sig.comb_cap,args_out_sig.ic_cap,args_out_sig.an_cap];

ABR_fs = 48828;
dur = 0.0310;

ihc_grades = [50,10,3,0.1,0.01]/100;

Fs = 100e3;
fs=Fs;
fc = 4000;
fm = 103;
modDepth = 1;
% dur = length(pitches(:,1))/fs; %use same dur
signallevel = 90;
Es=signallevel;
amp=20e-6*10.^(Es/20);
%dur = 1;

[sam_tone, t]= make_SAM(fc, fm, fs, modDepth, dur, amp,[],[]);
[r25, t] = make_RAM(fc, fm, fs, modDepth, 25, dur, amp, [], []);
[r50, t] = make_RAM(fc, fm, fs, modDepth, 50, dur, amp, [], []);

%should have all my stims
all_stims = [sam_tone,r50',r25'];
% all_stims = sam_tone;

for r = 1:size(all_stims,2)
    
    modelParams.cihc = ones(length(CF),1); %healthy
    fprintf('\n STIMULUS %i of %i',r,size(all_stims,2));
    input = all_stims(:,r);
    input = input';
    modelParams.dur = length(input)/fs;
    fprintf('\n    Normal Run');
    args_out_sig = make_EvokedPotential_new(input,CFs,unitary_resp,ABR_fs,modelParams,ihc_loss_flag); % add ABR_fs as o/p, CFs as i/p
    grand_envs_n(:,r) = args_out_sig;
    
    for g = 1:length(ihc_grades)
        fprintf('\n    Impaired Run %i of %i', g, length(ihc_grades));
        modelParams.cihc = ihc_grades(g).*ones(length(CF),1);
        args_out_sig = make_EvokedPotential_new(input,CFs,unitary_resp,ABR_fs,modelParams,ihc_loss_flag); % add ABR_fs as o/p, CFs as i
        grand_envs_i(g,:,r) =  args_out_sig;
    end
    
end

timestr = [num2str(yyyymmdd(datetime)),num2str(hour(datetime)),num2str(minute(datetime))];
cd(simdir);
save(['sim_',timestr,'.mat']);
cd(cwd);

