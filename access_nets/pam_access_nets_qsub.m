function ber = pam_access_nets_qsub(M, fiberLengthKm, wavelengthnm, ModBWGHz, amplified, RecBWGHz, Ntaps, ENOB, ros)
%% Calculate BER of IM-DD system using M-PAM
% - Equalization is done using a fractionally spaced linear equalizer
% Simulations include modulator, fiber, optical amplifier (optional) characterized 
% only by gain and noise figure, optical bandpass filter, antialiasing 
% filter, sampling, and linear equalization

addpath ../f % general functions
addpath ../apd % for PIN photodetectors
addpath ../apd/f

%% Transmit power swipe
% Tx. PtxdBm = -30:-20;
% PtxdBm = -25:-20; % transmitter power range
Tx.PtxdBm = -30:0.5:-10; % transmitter power range

%% Simulation parameters
sim.Rb = 25e9;    % bit rate in bits/sec
sim.Nsymb = 2^16; % Number of symbols in montecarlo simulation
sim.ros.txDSP = 1; % oversampling ratio transmitter DSP (must be integer). DAC samping rate is sim.ros.txDSP*mpam.Rs
% For DACless simulation must make Tx.dsp.ros = sim.Mct and DAC.resolution = Inf
sim.ros.rxDSP = str2double(ros); % oversampling ratio of receiver DSP. If equalization type is fixed time-domain equalizer, then ros = 1
sim.Mct = 2*5;      % Oversampling ratio to simulate continuous time. Must be integer multiple of sim.ros.txDSP and numerator of sim.ros.rxDSP
sim.BERtarget = 1.8e-4; 
sim.Ndiscard = 512; % number of symbols to be discarded from the begning and end of the sequence
sim.N = sim.Mct*sim.Nsymb; % number points in 'continuous-time' simulation
sim.Modulator = 'MZM'; % 'MZM' or 'DML'

sim.save = true;
 
%% Simulation control
sim.preAmp = amplified==1;        % optical amplifier
sim.WhiteningFilter = false;
sim.OptimizeGain = false;
sim.preemphasis = false; % preemphasis to compensate for transmitter bandwidth limitation
sim.preemphRange = 25e9; % preemphasis range
sim.mzm_predistortion = 'levels'; % predistortion to compensate MZM nonlinearity {'none': no predistortion, 'levels': only PAM levels are predistorted, 'analog': analog waveform is predistorted (DEPRECATED)}
sim.RIN = true; % include RIN noise in montecarlo simulation
sim.phase_noise = true; % whether to simulate laser phase noise
sim.PMD = false; % whether to simulate PMD
sim.quantiz = true; % whether quantization is included
sim.stopSimWhenBERReaches0 = true; % stop simulation when counted BER reaches 0

% Control what should be plotted
sim.Plots = containers.Map();
sim.Plots('BER') = 1;
sim.Plots('DAC output') = 0;
sim.Plots('Optical eye diagram') = 0;
sim.Plots('Received signal eye diagram') = 0;
sim.Plots('Signal after equalization') = 1;
sim.Plots('Equalizer') = 1;
sim.Plots('Electronic predistortion') = 0;
sim.Plots('Adaptation MSE') = 0;
sim.Plots('Channel frequency response') = 0;
sim.Plots('OSNR') = 0;
sim.Plots('Decision errors') = 0;
sim.Plots('Received signal optical spectrum') = 0;
sim.Plots('PAM levels MZM predistortion') = 0;
sim.shouldPlot = @(x) sim.Plots.isKey(x) && sim.Plots(x);

%% Pulse shape
Tx.pulse_shape = select_pulse_shape('rect', sim.ros.txDSP);
% pulse_shape = select_pulse_shape('rrc', sim.ros.txDSP, 0.5, 6);

%% M-PAM
% PAM(M, bit rate, leve spacing : {'equally-spaced', 'optimized'}, pulse
% shape: struct containing properties of pulse shape 
mpam = PAM(str2double(M), sim.Rb, 'equally-spaced', Tx.pulse_shape);

%% Time and frequency
sim.fs = mpam.Rs*sim.Mct;  % sampling frequency in 'continuous-time'
[sim.f, sim.t] = freq_time(sim.N, sim.fs);

%% ==========================  Transmitter ================================
%% DAC
Tx.DAC.fs = sim.ros.txDSP*mpam.Rs; % DAC sampling rate
Tx.DAC.ros = sim.ros.txDSP; % oversampling ratio of transmitter DSP
Tx.DAC.resolution = Inf; % DAC effective resolution in bits
Tx.DAC.rclip = 0;
Tx.DAC.filt = design_filter('bessel', 5, 0.7*mpam.Rs/(sim.fs/2)); % DAC analog frequency response
% juniperDACfit = [-0.0013 0.5846 0];
% Tx.DAC.filt.H = @(f) 1./(10.^(polyval(juniperDACfit, abs(f*sim.fs/1e9))/20)); % DAC analog frequency response

%% Laser
% Laser constructor: laser(lambda (nm), PdBm (dBm), RIN (dB/Hz), linewidth (Hz), frequency offset (Hz))
% lambda : wavelength (nm)
% PdBm : output power (dBm)
% RIN : relative intensity noise (dB/Hz)
% linewidth : laser linewidth (Hz)
% freqOffset : frequency offset with respect to wavelength (Hz)
wavelength = str2double(wavelengthnm)*1e-9;
Tx.Laser = laser(wavelength, 0, -150, 0.2e6, 0);
Tx.Laser.alpha = 0;

%% Modulator
Tx.Mod.rexdB = -15;  % extinction ratio in dB. Defined as Pmin/Pmax
Tx.Mod.type = sim.Modulator; 
Tx.Vgain = 1; % Gain of driving signal
Tx.VbiasAdj = 1; % adjusts modulator bias
Tx.Mod.Vbias = 0.5; % bias voltage normalized by Vpi
Tx.Mod.Vswing = 1;  % normalized voltage swing. 1 means that modulator is driven at full scale
Tx.Mod.BW = str2double(ModBWGHz)*1e9; % DAC frequency response includes DAC + Driver + Modulator
Tx.Mod.filt = design_filter('two-pole', Tx.Mod.BW, sim.fs);
Tx.Mod.H = Tx.Mod.filt.H(sim.f/sim.fs);
Tx.Mod.alpha = 0;

%% ============================= Fiber ====================================
% fiber(Length in m, anonymous function for attenuation versus wavelength
% (default: att(lamb) = 0 i.e., no attenuation), anonymous function for 
% dispersion versus wavelength (default: SMF28 with lamb0 = 1310nm, S0 = 0.092 s/(nm^2.km))
SMF = fiber(str2double(fiberLengthKm)*1e3); 
DCF = fiber(0, @(lamb) 0, @(lamb) -0.1*(lamb-1550e-9)*1e3 - 40e-6); 
DCF.L = SMF.L*SMF.D(wavelength)/DCF.D(wavelength);

Fibers = [SMF DCF];
% Fibers = SMF;

linkAttdB = SMF.att(Tx.Laser.wavelength)*SMF.L/1e3...
    + DCF.att(Tx.Laser.wavelength)*DCF.L/1e3;

%% ========================== Amplifier ===================================
% Constructor: OpticalAmplifier(Operation, param, Fn, Wavelength)
% - Opertation: either 'ConstantOutputPower' or 'ConstantGain'
% - param: GaindB if Operation = 'ConstantGain', or outputPower
% if Operation = 'ConstantOutputPower'
% - Fn:  noise figure in dB
% - Wavelength: operationl wavelength in m
Rx.OptAmp = OpticalAmplifier('ConstantOutputPower', 0, 5, Tx.Laser.wavelength); 
% Rx.OptAmp = OpticalAmplifier('ConstantGain', 20, 5, Tx.Laser.wavelength);

%% ============================ Receiver ==================================
%% Photodiodes

% PIN
% pin(R, Id, BW)
if amplified == 2
    % APD
    % apd(GaindB, ka, BW, R, Id)
    Rx.PD = apd(12, 0.5, str2double(RecBWGHz)*1e9, 1, 10e-9);
else
    Rx.PD = pin(1, 10e-9, Inf);
end

%% TIA-AGC
% One-sided thermal noise PSD
Rx.N0 = (30e-12).^2; 

%% Receiver DSP
Rx.filtering = 'antialiasing'; % {'antialiasing' or 'matched'}

%% ADC for direct detection case
Rx.ADC.fs = sim.ros.rxDSP*mpam.Rs;
Rx.ADC.ros = sim.ros.rxDSP;
Rx.ADC.filt = design_filter('butter', 5, (Rx.ADC.fs/2)/(sim.fs/2)); % Antialiasing filter
Rx.ADC.ENOB = str2double(ENOB); % effective number of bits. Quantization is only included if quantiz = true and ENOB ~= Inf
Rx.ADC.rclip = 0;

%% Equalizer
% Terminology: TD = time domain, SR = symbol-rate, LE = linear equalizer
Rx.eq.ros = sim.ros.rxDSP;
Rx.eq.type = 'adaptive TD-LE'; %'fixed td-sr-le' or 'adaptive td-le'
Rx.eq.Ntaps = str2double(Ntaps); % number of taps
Rx.eq.mu = 1e-3; % adaptation ratio
Rx.eq.Ntrain = 1e4; % Number of symbols used in training (if Inf all symbols are used)
Rx.eq.Ndiscard = [1.2e4 128]; % symbols to be discard from begining and end of sequence due to adaptation, filter length, etc


%% optimize APD gain and PAM level spacing
if ~isa(Rx.PD,'pin') && isfield(sim, 'OptimizeGain') && sim.OptimizeGain
    [Rx.PD.Gain, mpam] = Rx.PD.optGain(mpam, Tx, SMF, Rx, sim);
    fprintf('Optimal APD Gain = %.2f (%2.f dB)\n', Rx.PD.Gain, Rx.PD.GaindB);
    % if mpam.level_spacing = 'optimized', then Apd.optGain returns mpam 
    % with optimal level spacing
elseif mpam.optimize_level_spacing  %% Level Spacing Optimization
    % Optimize levels using Gaussian approximation
    [~, mpam] = Rx.PD.optimize_PAM_levels(Rx.PD.Gain, mpam, Tx, SMF, Rx, sim);
    mpam = mpam.norm_levels();
end
%% Generate summary
% generate_summary(mpam, Tx, Fibers, Rx.OptAmp, Rx, sim);

% check if there are enough symbols to perform simulation
Ndiscarded = sum(Rx.eq.Ndiscard) + 2*sim.Ndiscard;
assert(Ndiscarded < sim.Nsymb, 'There arent enough symbols to perform simulation. Nsymb must be increased or Ndiscard must be reduced')
fprintf('%d (2^%.2f) symbols will be used to calculated the BER\n', sim.Nsymb - Ndiscarded, log2(sim.Nsymb - Ndiscarded));

%% Run simulation
%% Calculate BER
Ptx = dBm2Watt(Tx.PtxdBm); % Transmitted power

ber.gauss = zeros(size(Ptx));
ber.count = zeros(size(Ptx));
OSNRdB = zeros(size(Ptx));
for k = 1:length(Ptx)
    Tx.Ptx = Ptx(k);
% Tx.Ptx = Ptx(1);
% gains = 8:13;
% for k = 1:length(gains)  
%     Rx.PD.GaindB = gains(k);
    % Montecarlo simulation
    [ber.count(k), ber.gauss(k), OSNRdB(k), Rx] = ber_pam_montecarlo(mpam, Tx, Fibers, Rx, sim);
    % only equalizer and noiseBW fields are changed in struct Rx
    
    if sim.stopSimWhenBERReaches0 && ber.count(k) == 0
        break;
    end 
end

%% Plots
if sim.shouldPlot('BER') && length(ber.count) > 1
    figure(1), hold on, box on
    if sim.preAmp
        hline(1) = plot(OSNRdB, log10(ber.count), '-o', 'LineWidth', 2, 'DisplayName', 'Counted');
        hline(2) = plot(OSNRdB, log10(ber.gauss), '-', 'Color', get(hline(1), 'Color'), 'LineWidth', 2, 'DisplayName', 'Gaussian Approx.');
        hline(3) = plot(OSNRdB(OSNRdB ~= 0), log10(pam_ber_from_osnr(mpam.M, OSNRdB(OSNRdB ~= 0), Rx.noiseBW)), '--k', 'LineWidth', 2, 'DisplayName', 'Sig-spont & noise enhancement');
        hline(3) = plot(OSNRdB(OSNRdB ~= 0), log10(pam_ber_from_osnr(mpam.M, OSNRdB(OSNRdB ~= 0), mpam.Rs/2)), ':k', 'LineWidth', 2, 'DisplayName', 'Sig-spont limit');
        legend('-DynamicLegend')
        axis([OSNRdB(1) OSNRdB(find(OSNRdB ~= 0, 1, 'last')) -8 0])
        xlabel('OSNR (dB)', 'FontSize', 12)
    else
        PrxdBm = Tx.PtxdBm - linkAttdB;
%         PrxdBm = gains;
        hline(1) = plot(PrxdBm, log10(ber.count), '-o', 'LineWidth', 2, 'DisplayName', 'Counted');
        hline(2) = plot(PrxdBm, log10(ber.gauss), '-', 'Color', get(hline(1), 'Color'), 'LineWidth', 2, 'DisplayName', 'Gaussian Approximation');
        legend('-DynamicLegend')
        axis([PrxdBm(1) PrxdBm(end) -8 0])
        xlabel('Received power (dBm)', 'FontSize', 12)
%         xlabel('APD gain (dB)','FontSize',12)
    end
end
ylabel('log_{10}(BER)', 'FontSize', 12)
set(gca, 'FontSize', 12)

%% save results
if sim.save
    folder = sprintf('results/%dPAM/',mpam.M);
    if ~exist(folder,'dir')
        mkdir(folder)
    end
    try
        filename = sprintf('results/%dPAM/L=%skm/PAM_BER_lamb=%snm_ModBW=%dGHz_',...
            mpam.M,fiberLengthKm, wavelengthnm, Tx.Mod.BW*1e-9);
        if isa(Rx.PD,'pin')
            filename = [filename sprintf('amplified=%d_', sim.preAmp)];
        else % using apd
            filename = [filename 'amplified=apd_'];
        end
        if Rx.PD.BW < Inf
            filename = [filename sprintf('RecBW=%s_', RecBWGHz)];
        else
            filename = [filename 'RecBW=Inf_'];
        end
        filename = [filename sprintf('Ntaps=%d_ENOB=%d_ros=%.2f.mat', Rx.eq.Ntaps, Rx.ADC.ENOB, sim.ros.rxDSP)];

        filename = check_filename(filename)
        
    % delete large variables
    sim = rmfield(sim, 'f');
    sim = rmfield(sim, 't');
    Tx.Mod = rmfield(Tx.Mod, 'H');    
    save(filename)
    catch
        warning('error saving file, saving to temp.mat')
        if isfield(sim, 'f')
            sim = rmfield(sim, 'f');
            sim = rmfield(sim, 't');
            Tx.Mod = rmfield(Tx.Mod, 'H');  
        end
        save('temp.mat')
    end
end