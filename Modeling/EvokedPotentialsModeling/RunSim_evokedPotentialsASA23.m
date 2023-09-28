%Andrew Sivaprakasam
%Description: Main file to run to simulate evoked potentials using the
%BEZ2018 model. 

%Key Assumptions in this Version:
% -Using only CFs that are multiples of stimulus F0s  
% -EFR (in animal model) mostly dominated by periphery/AN responses
% -Single SR chosen,using mean of SRs collected in our lab

%tic
%% Clearing and Adding Paths
clear all, close all

addpath('Stimulus_Generation')
addpath('BEZ2018model/')
addpath('Functions')

cwd = pwd();
%where to save simulations
simdir = 'SimulatedData';

if ~isfolder(simdir)
    mkdir(simdir);
end

%% Compiling C Code
% 
% cd BEZ2018model
% mexANmodel
% cd ../ 
% 
% cd Functions
% mex gammatone_c.c
% cd ../

%% Model Parameter Initialization:

%modelParams.CF = 0;
F0 = 103;
CF = [1.5:0.5:13,16,20,4e3/F0,40,80,160,190]*F0; %make sure to simulate 4k chan
NFFT = 4000;
NW = 3;
dB_stim = 75;
spontrates = [55,40]; %[Normal, Impaired] to match in vivo mean

%Normal Model Params
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

%Impaired Params (just changing cohc/cihc):
% [cohc_impaired,cihc_impaired,~] = fitaudiogram2(CF,dB_loss,modelParams.species);
ihc_grades = [50,10,3,0.1,0.01]/100;

%% Stimuli Initialization:

load("pitch_harmonic_rank_alt.mat");
fc = 4000;
fm = 103;
fs = fs; %use what's defined in my stim matrix
modDepth = 1;
dur = length(pitches(:,1))/fs; %use same dur
amp = 1;
iters = 50;

[sam_tone, t]= make_SAM(fc, fm, fs, modDepth, dur, amp,[],[]);
[r25, t] = make_RAM(fc, fm, fs, modDepth, 25, dur, amp, [], []);
[r50, t] = make_RAM(fc, fm, fs, modDepth, 50, dur, amp, [], []);


%should have all my stims
all_stims = [pitches,sam_tone,r50',r25'];
% all_stims = sam_tone;

for r = 1:size(all_stims,2)
    
    modelParams.cihc = ones(length(CF),1); %healthy
    modelParams.spont = spontrates(1);
    fprintf('\n STIMULUS %i of %i',r,size(all_stims,2));
    input = all_stims(:,r);
    modelParams.dur = length(input)/fs;
    fprintf('\n    Normal Run');
    [psth_pos,psth_neg,~] = getAP_PSTH(input,fs,modelParams,CF);
    grand_envs_n(:,r) = sum(psth_neg,1)+sum(psth_pos,1);
    
    for g = 1:length(ihc_grades)
        fprintf('\n    Impaired Run %i of %i', g, length(ihc_grades));
        modelParams.cihc = ihc_grades(g).*ones(length(CF),1);
        modelParams.spont = spontrates(2);
        [psth_pos,psth_neg,psth_fs] = getAP_PSTH(input,fs,modelParams,CF);
        grand_envs_i(g,:,r) = sum(psth_neg,1)+sum(psth_pos,1);
    end
    
end

timestr = [num2str(yyyymmdd(datetime)),num2str(hour(datetime)),num2str(minute(datetime))];

cd(simdir);
save(['sim_',timestr,'.mat']);
cd(cwd);
%% quick filtering
% 
% band = [75,300];
% band = band/(psth_fs/2);
% 
% [b,a] = butter(6,band);
% 
% grand_env_alt_filt = filtfilt(b,a,grand_env_alt);
