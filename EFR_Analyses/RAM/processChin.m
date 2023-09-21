%Author: Andrew Sivaprakasam
%Updated: July, 2023
%Purpose: Script to import/plot/apply additional processing to RAM_EFR
%files (chin version)

%Helpful Info: Be sure to define datapath so Import data section works. 
%see my example setup_AS file.   


harmonics = 16;

fs = 8e3; %fs to resample to
t_win = [.2,.9]; %signal window, ignoring onset/offset effects
filts = [60,4000];
frames = round(t_win*fs);

%% Handles my Local File Structure/EXT drive

%% Import data
cwd = pwd;
cd(datapath)
datafile = {dir(fullfile(cd,search_str)).name};
load(datafile{1});
fname_out = [datafile{1}(1:end-4),'_matlab.mat'];
cd(cwd);
%% Data analysis & plotting:
fs_orig = data.Stimuli.RPsamprate_Hz;
fs = 8e3; %resample to 8kHz

if size(data.AD_Data.AD_All_V,2)==1 || size(data.AD_Data.AD_All_V,2)==2
    all_dat = cell2mat(data.AD_Data.AD_All_V{1}');
else
    all_dat = cell2mat(data.AD_Data.AD_All_V');
end
    all_dat = all_dat';

[b,a] = butter(4,filts./(fs_orig/2));
all_dat = filtfilt(b,a,all_dat);

all_dat = resample(all_dat,fs,round(fs_orig));
all_dat = all_dat(frames(1):frames(2),:);

pos = all_dat(:,1:2:end)*1e6/data.AD_Data.Gain; %+ polarity
neg = all_dat(:,2:2:end)*1e6/data.AD_Data.Gain; %- polarity

%% Get PLV spectra/Time domain waveform:

%params for random sampling with replacement
subset = 100;
k_iters = 30;

%only output things we want to look at
[f, ~, ~, PLV_env, ~, ~, T_env] = helper.getSpectAverage(pos,neg, fs, subset, k_iters);
t = (1:length(T_env))/fs;

%% Get Peaks

[PKS,LOCS] = helper.getPeaks(f,PLV_env,fmod,harmonics);

%% Plot:
blck = [0.25, 0.25, 0.25];
rd = [0.8500, 0.3250, 0.0980, 0.5];
figure;

%Spectral Domain
hold on;
title([subj,' | RAM - 25% Duty Cycle | ',condition],'FontSize',14);
plot(f,PLV_env,'Color',blck,'linewidth',1.5);
plot(LOCS,PKS,'*','Color',rd,'MarkerSize',10,'LineWidth',2);

hold off;
ylim([0,1])
ylabel('PLV','FontWeight','bold')
xlabel('Frequency(Hz)','FontWeight','bold')

%Time Domain
xstart = .6;
xend = .9;
ystart = 0.6;
yend = .9;

axes('Position',[xstart ystart xend-xstart yend-ystart])
box on
hold on
plot(t, T_env,'Color',blck, 'LineWidth',2);
xlim([0.3,.4]);
ylim([-2,2]);
yticks([-1,0,1])
xlabel('Time(s)','FontWeight','bold');
ylabel('Amplitude \muV','FontWeight','bold')
hold off

set(gcf,'Position',[1557 538 560 420])

%% Export:

suffix2 = [condition,'/',subj,'/Preprocessed'];

data_out = [prefix,suffix2];
if ~exist(data_out,'dir')
    mkdir(data_out);
end

cd(data_out);
fname = [subj,'_RAM_efr_chin_',condition];
print(gcf,[fname,'_figure'],'-dpng','-r300');
save(fname,'t','T_env','f','PLV_env','PKS','LOCS')
cd(cwd)