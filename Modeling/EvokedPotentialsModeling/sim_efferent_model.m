function [ihcout, hsrout, lsrout, icout, gain] = sim_efferent_model(x, cf, args)
% SIM_EFFERENT_MODEL(x, cf) simulates efferent-model response to row-vector
% stimulus x at CFs in vector cf. 
% 
% See detailed help via `help sim_efferent_model`, including information
% about input arguments, return values, and changes to the code over time
% (changelog).
% 
% Returned values are matrices of size (n_chan, n_sample), where 
% `n_sample=length(x)` and `n_chan=length(cf)`. The various output matrices 
% are described below:
%   1) Inner hair cell "voltage", in a.u.
%   2) High-spontaneous-rate auditory-nerve instantaneous rate, in sp/s
%   3) Low-spontaneous-rate auditory-nerve instantaneous rate, in sp/s
%   4) Inferior-colliculus (IC) band-enhanced rate, in sp/s
%   5) Time-varying cochlear gain factor (in [0, 1], where 0==no gain, 1==max gain)
%
% SIM_EFFERENT_MODEL(x, cf) passes evaluates the efferent model on the
% input sound-pressure waveform at particular CFs. 
%
% SIM_EFFERENT_MODEL(x, cf, species=2) runs the efferent model for a 
% species value of 2 (which corresponds to human tuning based on data from 
% Shera). Other model parameters, such as sampling rate or IC model
% parameters, are adjusted the same way by specifying a key-value
% combinations (e.g., ic_tau_e=0.5e-3 would set the inhibitory IC delay to
% 0.5 ms, or moc_cutoff=2.0 would set the MOC lowpass cutoff to 2 Hz). See
% below for more details about available parameters and their default
% values (which are always used unless otherwise specified).
%
% Changelog:
% Changes to the MATLAB/Mex model code are documented here, while changes
% to the model code itself are documented in a separate changelog file.
%
% - 08/11/2023, DRG
%	Cleaned up Mex files, updated documentation in this file, updated this 
%   file to include validator functions applied to inputs to ensure that
%   only inputs within a sensible range can be provided.
%
% - 07/27/2023, DRG
%   Added the fast power-law adaptation approximation based on a parallel
%   set of exponential adaptation processes. Enabled by passing
%   `powerlaw_mode=2`.
%
% - 07/10/2023, DRG
%   Corrected synthesis of fractional Gaussian noise, which was previously
%   using inappropriate values for noise variance from the 2018 model code,
%   instead of correct values from the 2014 model code.
%
% - 06/29/2023, DRG
%   Adjusted MOC lowpass filter cutoff value to 0.64 Hz, which should more
%   closely match the value used in the old single-channel efferent model.
%
% - 06/21/2023, DRG
%   Fixed a memory leak in the model due to incorrect freeing of
%   dynamically allocated memory in the Mex file.
%
% Arguments:
% - x: Vector containing input sound-pressure waveform (Pa)
%
% - cf: Vector containing characteristic frequencies for each channel 
%   in the simulation (Hz)
%
% - args.fs: Sampling rate of the simulation (Hz). Note that inputs must be 
%   sampled at this sampling rate. 
%
% - args.cohc: Outer-hair-cell "count/health" (in [0, 1], where 0==no 
%   contribution of OHCs, 1==maximum contribution of OHCs)
%
% - args.cihc: Inner-hair-cell "count/health" (in [0, 1], where 0==no IHC 
%   amplitude, 1==normal IHC amplitude)
%
% - args.species: Which species to simulate in the basilar membrane/inner 
%   hair cell stage, 1==cat, 2==human (Shera), 3== (Moore and Glasberg)
%
% - args.powerlaw_mode: Whether to use true power-law adaptation
%   (powerlaw_mode == 1) or an approximate power-law implementation using a
%   set of 100 parallel exponential adaptation processes with time
%   constants fit computationally to match true power-law adaptation
%   (powerlaw_mode == 2). 
%
% - args.dur_settle: How long to simulate responses to silence before
%   simulating a response to input time-pressure waveform (s). Default of 
%   0.01 s (10 ms). This duration of time gives the various dynamic stages
%   of the model (e.g., AN adaptation, efferent gain control) a chance to
%   "settle in" to a more steady-state response regime before simulating
%   the stimulus response. If this is set to too short an interval, you may
%   see some weird response features at simulation onset.
%
% - args.ic_tau_e: Excitatory time constant in IC stage (s)
%
% - args.ic_tau_i: Inhibitory time constant in IC stage (s)
%
% - args.ic_delay: Inhibitory delay time in IC stage (s)
%
% - args.ic_amp: Excitatory strength in IC stage
%
% - args.ic_inh: Inhibitory strength in IC stage 
%
% - args.moc_cutoff: Cutoff of the lowpass filter used in the MOC stage
%   (Hz). The default value of 0.64 Hz yields a filter that matches that
%   used in the older single-channel efferent model (i.e., it produces a
%   "decay constant" of exp(-2pi * 0.64/100e3) ~= 1-3.9998e-5, which 
%   matches the constant used in the "old efferent" code, see Farhadi et 
%   al. 2023).
%
% - args.moc_beta_wdr: "beta" parameter in the MOC input-output
%   nonlinearity for the wide-dynamic-range MOC pathway (a.u.)
%
% - args.moc_offset_wdr: "offset" parameter in the MOC input-output
%   nonlinearity for the wide-dynamic-range MOC pathway (a.u.)
%
% - args.moc_beta_ic: "beta" parameter in the MOC input-output
%   nonlinearity for the IC MOC pathway (a.u.)
%
% - args.moc_offset_ic: "offset" parameter in the MOC input-output
%   nonlinearity for the IC MOC pathway (a.u.)
%
% - args.moc_weight_wdr: Scalar value multiplied with lowpass-filtered 
%   wide-dynamic-range MOC pathway signal before signal is passed through 
%   MOC input-output nonlinearity
%
% - args.moc_weight_ic: Scalar value multiplied with lowpass-filtered IC
%   MOC pathway signal before signal is passed through MOC input-output 
%   nonlinearity
%
% - args.moc_width_wdr: "Width" of the wide-dynamic-range cross-channel
%   "spread" (octaves). For example, a value of one octave means that each
%   wide-dynamic-range MOC signal will "spread" to all channels that have
%   CFs that fall within a band centered on the CF with a width of one
%   octave (within +/- one-half octave)
%
% - args.noiseType: Integer value determining whether we use empty matrices
%   (noiseType == -1), matrices of "frozen" fractional Gaussian noise 
%   (noiseType == 0), or matrices of "fresh" fractional Gaussian noise based
%   on the current global RNG state (noiseType == 1) as inputs for the
%   noise governing the stochastic behavior of the power-law synapse in the
%   auditory-nerve model. 
arguments
    x (1,:)
    cf (1,:) {mustBeGreaterThanOrEqual(cf, 100.0), mustBeLessThanOrEqual(cf, 40e3)}
    args.fs {mustBeGreaterThanOrEqual(args.fs,50e3), mustBeLessThanOrEqual(args.fs,200e3)} = 100e3
    args.cohc {mustBeGreaterThanOrEqual(args.cohc, 0.0), mustBeLessThanOrEqual(args.cohc, 1.0)} = 1.0
    args.cihc {mustBeGreaterThanOrEqual(args.cihc, 0.0), mustBeLessThanOrEqual(args.cihc, 1.0)} = 1.0
    args.species {mustBeMember(args.species, [1, 2, 3])} = 2
	args.powerlaw_mode {mustBeMember(args.powerlaw_mode, [1, 2])} = 2
    args.ic_tau_e {mustBeGreaterThan(args.ic_tau_e, 0.0)} = 1e-3
    args.ic_tau_i {mustBeGreaterThan(args.ic_tau_i, 0.0)} = 2e-3
    args.ic_delay {mustBeGreaterThanOrEqual(args.ic_delay, 0.0)} = 1e-3
    args.ic_amp {mustBeGreaterThan(args.ic_amp, 0.0)} = 1.0
    args.ic_inh {mustBeGreaterThanOrEqual(args.ic_inh, 0.0)} = 1.0
    args.moc_cutoff {mustBeGreaterThanOrEqual(args.moc_cutoff, 0.0)} = 0.64
    args.moc_beta_wdr {mustBeGreaterThanOrEqual(args.moc_beta_wdr, 0.0)} = 0.01
    args.moc_offset_wdr {mustBeGreaterThanOrEqual(args.moc_offset_wdr, 0.0)} = 0.0 
    args.moc_beta_ic {mustBeGreaterThanOrEqual(args.moc_beta_ic, 0.0)} = 0.01
    args.moc_offset_ic {mustBeGreaterThanOrEqual(args.moc_offset_ic, 0.0)} = 0.0
    args.moc_weight_wdr {mustBeGreaterThanOrEqual(args.moc_weight_wdr, 0.0)} = 2.0
    args.moc_weight_ic {mustBeGreaterThanOrEqual(args.moc_weight_ic, 0.0)} = 8.0
    args.moc_width_wdr {mustBeGreaterThanOrEqual(args.moc_width_wdr, 0.0)} = 0
    args.noiseType {mustBeMember(args.noiseType, [-1, 0, 1])} = 1
end

% Determine number of channels and samples
n_chan = length(cf);
n_sample = length(x);

% Synthesize fractional Gaussian noise
if args.noiseType == -1     
	% Use matrix of zeros
    ffGn_lsr = zeros(n_chan, n_sample);
    ffGn_hsr = zeros(n_chan, n_sample);
else                    
	% Synthesize fGn noise
    ffGn_lsr = zeros(n_chan, n_sample);
    ffGn_hsr = zeros(n_chan, n_sample);
    for ii=1:n_chan
        ffGn_lsr(ii, :) = ffGn(n_sample, 1/args.fs, 0.9, args.noiseType, 0.1, 3.0);
        ffGn_hsr(ii, :) = ffGn(n_sample, 1/args.fs, 0.9, args.noiseType, 100.0, 200.0);
    end
end

% Call Mex wrapper for efferent model
[ihcout, hsrout, lsrout, icout, gain] = sim_efferent_model_mex( ...
    x, ...
    ffGn_hsr, ...
    ffGn_lsr, ...
    cf, ...
    n_chan, ...
    1/args.fs, ...
    args.cohc, ...
    args.cihc, ...
    args.species, ...
    args.ic_tau_e, ...
    args.ic_tau_i, ...
    args.ic_delay, ...
    args.ic_amp, ...
    args.ic_inh, ...
    args.moc_cutoff, ...
    args.moc_beta_wdr, ...
    args.moc_offset_wdr, ...
    args.moc_beta_ic, ...
    args.moc_offset_ic, ...
    args.moc_weight_wdr, ...
    args.moc_weight_ic, ...
    args.moc_width_wdr, ...
	args.powerlaw_mode ...
    );

end