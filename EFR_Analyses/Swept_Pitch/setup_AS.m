%Here's where you can define your own parameters for input/output
%directories.
addpath(pwd);

close all;
clear;

condition = 'Carboplatin';
protocol = 'EFR_sweptPitch';
subj = 'Q403_awake';
swp_rnk = '_SweptRank';
% swp_rnk = 'REV_SweptRank';
freq = '103';
noise = '_20_';

search_str = ['p*',swp_rnk,'*',freq,'*',noise,'*.mat'];

fmod = 103;

uname = 'sivaprakasaman';
prefix = ['/media/',uname,'/AndrewNVME/Pitch_Study/PilotModulationEnhancementARO2024/',protocol,'/Chin/'];
suffix = [condition,'/',subj];

datapath = [prefix,suffix];
cwd = pwd;
cd(datapath)
datafile = {dir(fullfile(cd,search_str)).name};

write_out = [subj,'_',condition,'_',swp_rnk,'_',freq,'_',noise];

processChin;