clear 
close all
import mlreportgen.dom.*
import mlreportgen.report.*

report_name = ['Cat modeling results',date,'.pdf'];
rpt = Document(report_name, 'pdf'); %initialize report document as pdf
open(rpt);
pm = rpt.CurrentPageLayout;
%set page header dimensions
pm.PageMargins.Top = '0.1in';
pm.PageMargins.Header = '0.1in';
pm.PageMargins.Bottom = '0.1in';
pm.PageMargins.Footer = '0.1in';
pm.PageMargins.Left = '0.2in';
pm.PageMargins.Right = '0.2in';
name='human modeling results';
signallevels=40:10:90;
SNRlevels=-10:10:30;

signallevels=[60];
SNRlevels = [10];

for i=1:length(signallevels)
    for j=1:length(SNRlevels)
% 
signallevel=signallevels(i);
SNR=SNRlevels(j);

    h = Heading(2, ['signal level=',num2str(signallevel),'  SNR= ',num2str(SNR)]);
    b = Border();
    b.BottomStyle = 'single';
    b.BottomColor = 'LightGray';
    b.BottomWidth = '1pt';
    h.Style = [h.Style {Color('Darkorange'), b}];
    append(rpt,h);
% x=0:0.1:1;
% y=x.^2;
images = {};
% img=plot(x,y);


% 
% imgtype = '-dsvg';
%         imgname = [name '.svg'];
%         print(imgtype, imgname);
%         img = Image(imgname);
%         delete(gcf) %delete plot figure window
%         images = [images {img}];



        
        

%%
% [ihcout, hsr, lsr, ic, gain] = sim_efferent_model(...
% 	x,...
% 	cf,...
% 	moc_weight_wdr=16.0,..
% 	moc_weight_ic=8.0,...
% 	moc_width_wdr=0.5,...
% 	powerlaw_mode=2 ...
% );
 

fc=503;
fm=103;
modDepth=1;
% duty=25;
dur=1;
Es=signallevel;
fs=1e5;
amp=20e-6*10.^(Es/20);

signal=make_SAM(fc, fm, fs, modDepth, dur, amp) ;
% signalrms=20e-6*10.^(Es/20);
% signal=(signalrms/(rms(signal))).*signal;

noise40=ltass_noise(fs,signallevel-SNR,1,fs);
% noise60=ltass_noise(fs,60,1,fs);

SAM_503_0dB=signal';
SAM_503_40dB=signal'+noise40;
% SAM_503_60dB=signal'+noise60;

stiminput=SAM_503_40dB;
n_cf=50;%
CFs=logspace(log10(101.0), log10(2000.0), n_cf);
hsr_0dB_pop =zeros(n_cf,fs);
hsr_40dB_pop =zeros(n_cf,fs);

lsr_0dB_pop =zeros(n_cf,fs);
lsr_40dB_pop =zeros(n_cf,fs);

gain_0dB_pop =zeros(n_cf,fs);
gain_40dB_pop =zeros(n_cf,fs);

ic_0dB_pop =zeros(n_cf,fs);
ic_40dB_pop =zeros(n_cf,fs);

% cell(20,n_cf,fs);
% ihcout_40dB_pop = cell(20,n_cf,fs);
reps=2;

fs_efr = 8e3;

%AN and IC weights, make sure sum to 1
wts = [.4,.6];

parfor r=1:reps
    [p_ihcout_s, p_hsr_s, p_lsr_s, p_ic_s, p_gain_s] = sim_efferent_model(SAM_503_0dB,CFs);
    [n_ihcout_s, n_hsr_s, n_lsr_s, n_ic_s, n_gain_s] = sim_efferent_model(-SAM_503_0dB,CFs);
    
    %reshape to a new sample rate, only considering HSRs for now. 
    ap_psth_env = (p_hsr_s+n_hsr_s)/2;
    ap_psth_env = resample(ap_psth_env',fs_efr,fs);
    
    ap_psth_env_ic = (p_ic_s+n_ic_s)/2;
    ap_psth_env_ic = resample(ap_psth_env_ic',fs_efr,fs);

    %normalized
    an_cap = sum(ap_psth_env,2);
    an_cap = an_cap/max(an_cap);
    
    ic_cap = sum(ap_psth_env_ic,2);
    ic_cap = ic_cap/max(ic_cap);

    %TODO: convolve with click later???
    comb_cap = wts(1).*an_cap + wts(2).*ic_cap;

    [ihcout_40dB, hsr_40dB, lsr_40dB, ic_40dB, gain_40dB] = sim_efferent_model(SAM_503_40dB,CFs);

    ic_0dB_pop=ic_0dB_pop+ ic_0dB;
    ic_40dB_pop =ic_40dB_pop+ ic_40dB;
    
    hsr_0dB_pop =hsr_0dB_pop+hsr_0dB;
    hsr_40dB_pop =hsr_40dB_pop+hsr_40dB;
    
    lsr_0dB_pop =lsr_0dB_pop+lsr_0dB;
    lsr_40dB_pop =lsr_40dB_pop+lsr_40dB;
    
    gain_0dB_pop =gain_0dB_pop+gain_0dB;
    gain_40dB_pop =gain_40dB_pop+gain_40dB;

end

%%
%% 
     figure('Renderer', 'painters', 'Position', [10 10 800 500])
 

 hsr_0dB_average_reps=(hsr_0dB_pop/reps);
%  hsr_0dB_average_reps=squeeze(hsr_0dB_average_reps);

 hsr_40dB_average_reps=(hsr_40dB_pop/reps);
%  hsr_40dB_average_reps=squeeze(hsr_40dB_average_reps);


 lsr_0dB_average_reps=(lsr_0dB_pop/reps);
%  lsr_0dB_average_reps=squeeze(lsr_0dB_average_reps);

 lsr_40dB_average_reps=(lsr_40dB_pop/reps);
%  lsr_40dB_average_reps=squeeze(lsr_40dB_average_reps);


 gain_0dB_average_reps=(gain_0dB_pop/reps);
%  gain_0dB_average_reps=squeeze(gain_0dB_average_reps);

 gain_40dB_average_reps=(gain_40dB_pop/reps);
%  %  gain_40dB_average_reps=squeeze(gain_40dB_average_reps);

  ic_0dB_average_reps=(ic_0dB_pop/reps);
%  ic_0dB_average_reps=squeeze(ic_0dB_average_reps);

 ic_40dB_average_reps=(ic_40dB_pop/reps);
%  ic_40dB_average_reps=squeeze(ic_40dB_average_reps);



% hold on
ICaverage_0dB=mean(ic_0dB_average_reps(:,0.2*fs:0.9*fs),2);
HSraverage_0dB=mean(hsr_0dB_average_reps(:,0.2*fs:0.9*fs),2);
LSRavereg_0dB=mean(lsr_0dB_average_reps(:,0.2*fs:0.9*fs),2);
gainaverage_0dB=mean(gain_0dB_average_reps(:,0.2*fs:0.9*fs),2);

ICaverage_40dB=mean(ic_40dB_average_reps(:,0.2*fs:0.9*fs),2);
HSraverage_40dB=mean(hsr_40dB_average_reps(:,0.2*fs:0.9*fs),2);
LSRavereg_40dB=mean(lsr_40dB_average_reps(:,0.2*fs:0.9*fs),2);
gainaverage_40dB=mean(gain_40dB_average_reps(:,0.2*fs:0.9*fs),2);



subplot(1,4,1)
semilogx(CFs,HSraverage_0dB','LineWidth',3);
hold on
semilogx(CFs,HSraverage_40dB','LineWidth',3);
% legend('in silence', 'in noise','Location','northwestoutside')
title('HSR')
xlabel('CF (Hz)')
ylabel('Rate(spikes/sec)')
fontsize = 12;
set(gca,'fontsize',fontsize);
ylim([0,300])


subplot(1,4,2)
semilogx(CFs,gainaverage_0dB,'LineWidth',3);
hold on
semilogx(CFs,gainaverage_40dB,'LineWidth',3);
% legend('in silence', 'in noise')
title('Gain')
xlabel('CF (Hz)')
ylabel('Gain Factor')
fontsize = 12;
set(gca,'fontsize',fontsize);
ylim([0,1])


subplot(1,4,3)
semilogx(CFs,LSRavereg_0dB,'LineWidth',3);
hold on 
semilogx(CFs,LSRavereg_40dB,'LineWidth',3);
% legend('in silence', 'in noise')
title('LSR')
xlabel('CF (Hz)')
ylabel('Rate(spikes/sec)')
fontsize = 12;
set(gca,'fontsize',fontsize);
ylim([0,200])

subplot(1,4,4)
semilogx(CFs,ICaverage_0dB,'LineWidth',3);
hold on 
semilogx(CFs,ICaverage_40dB,'LineWidth',3);
% legend('in silence', 'in noise')
title('IC')
  xlabel('CF (Hz)')
ylabel('Rate(spikes/sec)')
fontsize = 12;
set(gca,'fontsize',fontsize);
ylim([0,100])
%%

     name = ['signal level=',num2str(signallevel),'  SNR= ',num2str(SNR)];
        imgtype = '-dsvg';
        imgname = [name '.svg'];
        print(imgtype, imgname);
        img = Image(imgname);
        delete(gcf) %delete plot figure window
        images = [images {img}];
        append(rpt, img);
        end
        end
close(rpt);