%Here's where you can define your own parameters for input/output
%directories.

close all;
clear;

condition = 'Baseline';
subj = 'Q427';
fmod = 103;

uname = 'sivaprakasaman';
prefix = ['/media/',uname,'/AndrewNVME/Pitch_Study/Pitch_Diagnostics_SH_AS/EFR_RAM/Chin/'];
suffix = [condition,'/',subj];

datapath = [prefix,suffix];

processChin;