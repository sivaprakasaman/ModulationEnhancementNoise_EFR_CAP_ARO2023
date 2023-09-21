%Here's where you can define your own parameters for input/output
%directories.

close all;
clear;

condition = 'Carboplatin';
subj = 'Q403';
search_str = 'p*RAM*503_60*.mat';
% search_str = 'p*FFR_REV_Swept*dbdrop__*.mat';
fmod = 103;

uname = 'sivaprakasaman';
prefix = ['/media/',uname,'/AndrewNVME/Pitch_Study/PilotModulationEnhancementARO2024/EFR_sweptPitch/Chin/'];
suffix = [condition,'/',subj];

datapath = [prefix,suffix];
processChin;