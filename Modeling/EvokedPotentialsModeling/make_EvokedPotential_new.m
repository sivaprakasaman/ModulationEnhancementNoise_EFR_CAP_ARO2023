function [hsr_conv_out] = make_EvokedPotential_new(stim,CFs,UR,UR_fs,modelParams,ihc_loss_flag)


% add 2018 and a switch to use either models
% [p_ihcout_s, p_hsr_s, p_lsr_s, p_ic_s, p_gain_s] = sim_efferent_model(stim,CFs);
% [n_ihcout_s, n_hsr_s, n_lsr_s, n_ic_s, n_gain_s] = sim_efferent_model(-stim,CFs);
% 

[p_hsr_s,n_hsr_s, psth_fs] = getAP_PSTH_cihc(stim,UR_fs,modelParams,CFs, ihc_loss_flag);
pp_hsr_s = sum(p_hsr_s,1);
np_hsr_s = sum(n_hsr_s,1);
sum_np_hsr = pp_hsr_s+np_hsr_s;




uresp=UR.f-mean(UR.f);
sum_np_hsr = resample(sum_np_hsr,UR_fs,psth_fs); %upsampling

%  psth_sum_sampled = sum_np_hsr(1,1:length(uresp));
hsr_conv_out=conv(sum_np_hsr, uresp);


end