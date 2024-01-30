function [psth_pos,psth_neg,psth_fs] = getAP_PSTH_cihc(input,input_fs,modelParams,CF, ihc_loss_flag)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

    spont = modelParams.spont;
    tabs = modelParams.tabs;
    trel = modelParams.trel;
    cohc = modelParams.cohc; %healthy
    cihc = modelParams.cihc; %healthy
    species = modelParams.species; % 1 for cat (2 for human with Shera et al. tuning; 3 for human with Glasberg & Moore tuning) %read up on this tuning
    noiseType = modelParams.noiseType; % 1 for variable fGn; 0 for fixed (frozen) fGn
    implnt = modelParams.implnt; % "0" for approximate or "1" for actual implementation of the power-law functions in the Synapse
    %stimdb = modelParams.stimdb;
    T = modelParams.dur; % stimulus duration in seconds
    nrep = modelParams.reps;
    Fs = modelParams.Fs;
    psthbinwidth = modelParams.psthbinwidth;
    buffer = modelParams.buffer;
    
    psthbins = round(psthbinwidth*Fs);
    %input = helper.gen_rescale(input, stimdb); %rescaling done outside
    pin = resample(input, Fs, input_fs); %upsampling
    pin = pin(1:fix(T*Fs));
    dt=1/modelParams.Fs;

    psth_pos = zeros(length(CF),buffer*(1/psthbinwidth)*ceil(T));
    psth_neg = zeros(length(CF),buffer*(1/psthbinwidth)*ceil(T));
    dc = 0.0006;
    
  parfor i = 1:length(CF)
        cf = CF(i);
        
        %ihc damage
        if ihc_loss_flag
            vihc_pos = model_IHC_BEZ2018_n(pin,cf,nrep,dt,buffer*ceil(T),cohc(i),cihc(i),species,dc);
            vihc_neg = model_IHC_BEZ2018_n(-pin,cf,nrep,dt,buffer*ceil(T),cohc(i),cihc(i),species,dc);

            [psthr_pos, ~, ~, ~, ~,~] = model_Synapse_BEZ2018_n(vihc_pos,cf,nrep,dt,noiseType,implnt,spont,tabs,trel);
            [psthr_neg, ~, ~, ~, ~,~] = model_Synapse_BEZ2018_n(vihc_neg,cf,nrep,dt,noiseType,implnt,spont,tabs,trel);
        else
            vihc_pos = model_IHC_BEZ2018(pin,cf,nrep,dt,buffer*ceil(T),cohc(i),cihc(i),species);
            vihc_neg = model_IHC_BEZ2018(-pin,cf,nrep,dt,buffer*ceil(T),cohc(i),cihc(i),species);

            [psthr_pos, ~, ~, ~, ~,~] = model_Synapse_BEZ2018(vihc_pos,cf,nrep,dt,noiseType,implnt,spont,tabs,trel);
            [psthr_neg, ~, ~, ~, ~,~] = model_Synapse_BEZ2018(vihc_neg,cf,nrep,dt,noiseType,implnt,spont,tabs,trel);
        end

        psth_pos(i,:) = sum(reshape(psthr_pos,psthbins,length(psthr_pos)/psthbins)); 
        psth_neg(i,:) = sum(reshape(psthr_neg,psthbins,length(psthr_neg)/psthbins));
        
   end
 
    psth_fs = Fs/psthbins;

end

