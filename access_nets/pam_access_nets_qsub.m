function [ber, sim, mpam, Tx, Fibers, Rx] = ...
pam_access_nets_qsub(RsGbd, M, fiberLengthKm, wavelengthnm, ModBWGHz, modulator, ...
levels, amplified, gain, RecBWGHz, equalization, BERtarget, PtxdBm, fc, savedata)
%% Calculate BER of IM-DD system using M-PAM
% - Equalization is done using a fractionally spaced linear equalizer
% Simulations include modulator, fiber, optical amplifier (optional) characterized 
% only by gain and noise figure, optical bandpass filter, antialiasing 
% filter, sampling, and linear equalization
% savedata - 0: no, 1: interactive, 2: autoskip, 3: rewrite

%% Transmit power swipe
if exist('PtxdBm','var')
    Tx.PtxdBm = PtxdBm;
else
    switch(M)
        case 2
            Tx.PtxdBm = -35:-25; % transmitter power range
    %         Tx.PtxdBm = -35:-20;
        case 4
            Tx.PtxdBm = -30:-15;
        case 8
    %         Tx.PtxdBm = -15:-5;
            Tx.PtxdBm = -25:-5;

        case 3
            Tx.PtxdBm = -30:-20;
    end
    % Tx.PtxdBm = -15;
end
Tx.PtxdBm = reshape(Tx.PtxdBm,[],1);
% Tx.PtxdBm = -22.5;
%% Simulation parameters
if M~=3
    sim.Rb = RsGbd*1e9*log2(M);    % bit rate in bits/sec
else
    sim.Rb = RsGbd*1e9*1.5;
end
sim.Nsymb = 2^17; % Number of symbols in montecarlo simulation
% For DACless simulation must make Tx.dsp.ros = sim.Mct and DAC.resolution = Inf
sim.ros.rxDSP = 1; % oversampling ratio of receiver DSP. If equalization type is fixed time-domain equalizer, then ros = 1
sim.Mct = 24;      % Oversampling ratio to simulate continuous time. Must be integer multiple of sim.ros.txDSP and numerator of sim.ros.rxDSP
sim.BERtarget = BERtarget;
sim.Ndiscard = 512; % number of symbols to be discarded from the begning and end of the sequence
sim.N = sim.Mct*sim.Nsymb; % number points in 'continuous-time' simulation
sim.Modulator = modulator; % 'MZM' or 'DML' or 'EAM'
if M~=3
    sim.L = 6; % de Bruijin sub-sequence length (ISI symbol length)
else
    sim.L = 4;
end
if ~isempty(regexp(equalization, 'precomp', 'once'))
    sim.ros.txDSP = 2; % oversampling ratio transmitter DSP (must be integer). DAC samping rate is sim.ros.txDSP*mpam.Rs
    sim.preemphasis = true;
    sim.preemphRange = 27e9;
    if strcmp(sim.Modulator,'MZM')
        sim.mzm_predistortion = 'analog'; % predistortion to compensate MZM nonlinearity {'none': no predistortion, 'levels': only PAM levels are predistorted, 'analog': analog waveform is predistorted (DEPRECATED)}
    else
        sim.mzm_predistortion = 'none';
    end
    sim.CDeq = true;
    numbers = regexp(equalization,'[0-9]');
    if ~isempty(numbers)
        sim.CDdist = str2num(equalization(numbers))*1e3;
    else
        sim.CDdist = [];
    end
        
    posteq = 'none';
else
    sim.preemphasis = false;
    sim.ros.txDSP = 1;
    posteq = equalization;
    if strcmp(sim.Modulator,'MZM')
        sim.mzm_predistortion = 'levels'; % predistortion to compensate MZM nonlinearity {'none': no predistortion, 'levels': only PAM levels are predistorted, 'analog': analog waveform is predistorted (DEPRECATED)}
    else
        sim.mzm_predistortion = 'none';
    end
end
sim.posteq = ~strcmp(posteq,'none');

if exist('savedata','var')
    sim.save = savedata;
else
    sim.save = false;
end
 
if sim.save
    params.lamb = wavelengthnm;
    
    params.RsGbd = RsGbd;
    params.ModBWGHz = ModBWGHz;
    params.RecBWGHz = RecBWGHz;
    params.Lkm = fiberLengthKm;
    params.modulator = modulator;
    amp = {'pin','soa','apd'};
    params.amp = amp{amplified+1};
    params.levels = levels;
    params.equalization = equalization;
    if exist('fc','var')
    params.fc = fc;
    end
    folder = 'C:/Users/Elaine/Documents/MATLAB/optical-comm/access_nets/results/';
    subfolder = sprintf('%dGbd/%dPAM/',RsGbd,M);
    if ~exist([folder subfolder],'dir')
        mkdir([folder subfolder])
    end
    filename = name(params,M,gain,BERtarget);
    if exist([folder filename],'file')
        switch sim.save
            case 1
                cont = input(sprintf('simulation data file %s already exists, rerun? (y/r/[n]): ',filename),'s');
                switch cont
                    case 'y'  % new file name
                        filename = check_filename([folder filename]);
                        ret = false;
                    case 'r' % rewrite old file
                        filename = [folder filename];
                        ret = false;
                    otherwise
                        ret = true;
                end
            case 2
                fprintf('simulation data file %s already exists, skipping\n',filename);
                ret = true;
            case 3
                fprintf('rewriting data file %s\n', filename);
                filename = [folder filename];
                ret = false;
        end
    else
        ret = false;
        filename = [folder filename];
    end
    if ret
        load([folder filename],'ber')
        return
    end
    disp(filename(length(folder)+1:end))
end

%% Simulation control
sim.preAmp = amplified==1;        % optical amplifier
sim.Apd = amplified==2;
sim.WhiteningFilter = ~sim.preemphasis;
sim.OptimizeGain = false;
sim.RIN = true; % include RIN noise in montecarlo simulation
sim.phase_noise = true; % whether to simulate laser phase noise
sim.PMD = false; % whether to simulate PMD
sim.quantiz = true; % whether quantization is included
sim.terminateWhenBERReaches0 = true; % stop simulation when counted BER reaches 0

% Control what should be plotted
sim.Plots = containers.Map();
sim.Plots('BER') = 1;
sim.Plots('DAC output') = 0;
if ~sim.posteq
    sim.Plots('Optical eye diagram') = 1;
    sim.Plots('Received signal eye diagram') = 1;
end
sim.Plots('Signal after equalization') = 0;
sim.Plots('Equalizer') = 0;
sim.Plots('Adaptation MSE') = 0;
sim.Plots('Channel frequency response') = 0;
sim.Plots('OSNR') = 0;
sim.Plots('Decision errors') = 0;
sim.Plots('Received signal optical spectrum') = 0;
sim.Plots('PAM levels MZM predistortion') = 0;
sim.Plots('Conditional PDF') = 1;
sim.shouldPlot = @(x) sim.Plots.isKey(x) && sim.Plots(x);

%% Pulse shape
Tx.pulse_shape = select_pulse_shape('rect', sim.ros.txDSP);
% pulse_shape = select_pulse_shape('rrc', sim.ros.txDSP, 0.5, 6);

%% M-PAM
% PAM(M, bit rate, leve spacing : {'equally-spaced', 'optimized'}, pulse
% shape: struct containing properties of pulse shape 

mpam = PAM(M, sim.Rb, levels, Tx.pulse_shape);

%% Time and frequency
sim.fs = mpam.Rs*sim.Mct;  % sampling frequency in 'continuous-time'
[sim.f, sim.t] = freq_time(sim.N, sim.fs);

%% ==========================  Transmitter ================================
%% DAC
Tx.DAC.fs = sim.ros.txDSP*mpam.Rs; % DAC sampling rate
Tx.DAC.ros = sim.ros.txDSP; % oversampling ratio of transmitter DSP
Tx.DAC.offset = sim.Mct*(1-1/sim.ros.txDSP)/2;
if sim.preemphasis
    Tx.DAC.resolution = 5; % DAC effective resolution in bits
    Tx.DAC.rclip = 0;
else
    Tx.DAC.resolution = Inf;
    Tx.DAC.rclip = 0;
end
Tx.DAC.filt = design_filter('bessel', 5, 0.7*Tx.DAC.fs/(sim.fs/2)); % DAC analog frequency response
% juniperDACfit = [-0.0013 0.5846 0];
% Tx.DAC.filt.H = @(f) 1./(10.^(polyval(juniperDACfit, abs(f*sim.fs/1e9))/20)); % DAC analog frequency response

%% Laser
% Laser constructor: laser(lambda (nm), PdBm (dBm), RIN (dB/Hz), linewidth (Hz), frequency offset (Hz))
% lambda : wavelength (nm)
% PdBm : output power (dBm)
% RIN : relative intensity noise (dB/Hz)
% linewidth : laser linewidth (Hz)
% freqOffset : frequency offset with respect to wavelength (Hz)
wavelength = wavelengthnm*1e-9;
% Tx.Laser = laser(wavelength, 0, -150, 0.2e6, 0);
Tx.Laser = laser(wavelength, 0);
Tx.Laser.alpha = 0;

%% Modulator
% Tx.Mod.rexdB = -10;  % extinction ratio in dB. Defined as Pmin/Pmax
Tx.Mod.rexdB = -inf;  % extinction ratio in dB. Defined as Pmin/Pmax
Tx.Mod.type = sim.Modulator; 
Tx.Vgain = 1; % Gain of driving signal
Tx.VbiasAdj = 1; % adjusts modulator bias
Tx.Mod.Vbias = 0.5; % bias voltage normalized by Vpi
Tx.Mod.Vswing = 1;  % normalized voltage swing. 1 means that modulator is driven at full scale
Tx.Mod.BW = ModBWGHz*1e9; % DAC frequency response includes DAC + Driver + Modulator
Tx.Mod.filt = design_filter('two-pole', Tx.Mod.BW, sim.fs);
Tx.Mod.H = Tx.Mod.filt.H(sim.f/sim.fs);
Tx.Mod.alpha = 0;

%% ============================= Fiber ====================================
% fiber(Length in m, anonymous function for attenuation versus wavelength
% (default: att(lamb) = 0 i.e., no attenuation), anonymous function for 
% dispersion versus wavelength (default: SMF28 with lamb0 = 1310nm, S0 = 0.092 s/(nm^2.km))
fiberlen = fiberLengthKm*1e3;
SMF = fiber(fiberlen); 
% DCF = fiber(0, @(lamb) 0, @(lamb) -0.1*(lamb-1550e-9)*1e3 - 40e-6); 
% DCF.L = SMF.L*SMF.D(wavelength)/DCF.D(wavelength);

% Fibers = [SMF DCF];
Fibers = SMF;

linkAttdB = SMF.att(Tx.Laser.wavelength)*SMF.L/1e3;%...
%     + DCF.att(Tx.Laser.wavelength)*DCF.L/1e3;

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

if sim.Apd
    % APD
    % apd(GaindB, ka, BW, R, Id)
    Rx.PD = apd(gain, 0.15, [RecBWGHz*1e9 300e9] , 0.7, 10e-9);
%     Rx.PD = apd(gain, 0.15, RecBWGHz*1e9 , 0.7, 10e-9);
else
    % PIN3
    % pin(R, Id, BW)
    Rx.PD = pin(1, 10e-9, Inf);
end

%% TIA-AGC
% One-sided thermal noise PSD
Rx.N0 = (30e-12).^2; 

%% Receiver DSP
if sim.posteq
    Rx.filtering = 'matched'; % {'antialiasing' or 'matched'}
else
    Rx.filtering = 'antialiasing';
end

%% ADC for direct detection case
if sim.posteq
    Rx.ADC.fs = sim.ros.rxDSP*mpam.Rs;
    Rx.ADC.ros = sim.ros.rxDSP;
    Rx.ADC.filt = design_filter('butter', 5, (Rx.ADC.fs/2)/(sim.fs/2)); % Antialiasing filter
    Rx.ADC.ENOB = 5; % effective number of bits. Quantization is only included if quantiz = true and ENOB ~= Inf
    Rx.ADC.rclip = 0;
else
    assert(sim.ros.rxDSP==1,'no oversampling without equalization')
    Rx.ADC.fs = mpam.Rs;
    Rx.ADC.ros = 1;
    % increase filter cutoff from antialiasing filter
    Rx.ADC.fc = fc;
    Rx.ADC.filt = design_filter('bessel', 5, fc*Rx.ADC.fs/(sim.fs/2));
    Rx.ADC.ENOB = Inf;
    Rx.ADC.rclip = 0;
end

%% Equalizer
% Terminology: TD = time domain, SR = symbol-rate, LE = linear equalizer
Ndiscarded = 2*sim.Ndiscard;
if strcmp(posteq,'none')
    Rx.eq.ros = 1;
    Rx.eq.type = 'none';
    Rx.eq.Ndiscard = [0 0];
else
    Rx.eq.ros = sim.ros.rxDSP;
    Rx.eq.type = posteq; %'fixed td-sr-le' or 'adaptive td-le' or 'none'
    Rx.eq.Ntaps = 15; % number of taps
    Rx.eq.mu = 1e-3; % adaptation ratio
    Rx.eq.Ntrain = 2e4; % Number of symbols used in training (if Inf all symbols are used)
    Rx.eq.Ndiscard = [2.2e4 128]; % symbols to be discard from begining and end of sequence due to adaptation, filter length, etc
    Ndiscarded = Ndiscarded + sum(Rx.eq.Ndiscard);
end


%% Generate summary
% generate_summary(mpam, Tx, Fibers, Rx.OptAmp, Rx, sim);

% check if there are enough symbols to perform simulation
assert(Ndiscarded < sim.Nsymb, 'There arent enough symbols to perform simulation. Nsymb must be increased or Ndiscard must be reduced')
% fprintf('%d (2^%.2f) symbols will be used to calculated the BER\n', sim.Nsymb - Ndiscarded, log2(sim.Nsymb - Ndiscarded));

%% Calculate BER
[ber, mpam, Rx.PD] = apd_ber(mpam, Tx, SMF, Rx.PD, Rx, sim);
% display(Rx.PD.GaindB)
% mpam.level_spacing = 'optimized';
% [ber_apd_eq, mpam, Rx.PD] = apd_ber(mpam, Tx, SMF, Rx.PD, Rx, sim);
% display(Rx.PD.GaindB)
%% Run simulation
% Ptx = dBm2Watt(Tx.PtxdBm); % Transmitted power
% 
% ber.count2 = zeros(size(Ptx));
% ber.gauss = zeros(size(Ptx));
% OSNRdB = zeros(size(Ptx));
% for k = 1:length(Ptx)
%     Tx.Ptx = Ptx(k);
%     % Montecarlo simulation
%     [ber.count2(k), ber.gauss(k), OSNRdB(k), Rx] = ber_pam_montecarlo(mpam, Tx, Fibers, Rx, sim);
%     % only equalizer and noiseBW fields are changed in struct Rx
%     
%     if sim.stopSimWhenBERReaches0 && ber.count2(k) == 0
%         break;
%     end 
% end


%% Plots
if sim.preAmp
    PrxdBm = OSNRdbB;
else
    PrxdBm = Tx.PtxdBm - linkAttdB;
end
ber.PrxdBm = PrxdBm;
% if sim.shouldPlot('BER') && length(ber.count) > 1
%     figure(1), hold on, box on
%     if sim.preAmp
%         hline(1) = plot(OSNRdB, log10(ber.count2), '-o', 'LineWidth', 2, 'DisplayName', 'Counted');
%         hline(2) = plot(OSNRdB, log10(ber.gauss), '-', 'Color', get(hline(1), 'Color'), 'LineWidth', 2, 'DisplayName', 'Gaussian Approx.');
%         hline(3) = plot(OSNRdB(OSNRdB ~= 0), log10(pam_ber_from_osnr(mpam.M, OSNRdB(OSNRdB ~= 0), Rx.noiseBW)), '--k', 'LineWidth', 2, 'DisplayName', 'Sig-spont & noise enhancement');
%         hline(3) = plot(OSNRdB(OSNRdB ~= 0), log10(pam_ber_from_osnr(mpam.M, OSNRdB(OSNRdB ~= 0), mpam.Rs/2)), ':k', 'LineWidth', 2, 'DisplayName', 'Sig-spont limit');
%         legend('-DynamicLegend')
%         axis([OSNRdB(1) OSNRdB(find(OSNRdB ~= 0, 1, 'last')) -8 0])
%         xlabel('OSNR (dB)', 'FontSize', 12)
%     else
%         PrxdBm = Tx.PtxdBm - linkAttdB;
% %         PrxdBm = gains;
%         hline(1) = plot(PrxdBm, log10(ber.count2), '-o', 'LineWidth', 2, 'DisplayName', 'Counted2');
%         hline(2) = plot(PrxdBm, log10(ber.gauss), '-', 'Color', get(hline(1), 'Color'), 'LineWidth', 2, 'DisplayName', 'Gaussian Approximation');
%         hline(3) = plot(PrxdBm, log10(ber.count), '-o', 'DisplayName', 'Counted');
%         hline(4) = plot(PrxdBm, log10(ber.awgn), '--', 'DisplayName', 'awgn');
%         hline(5) = plot(PrxdBm, log10(ber.enum), '--', 'DisplayName', 'enum');
%         legend('-DynamicLegend')
%         axis([PrxdBm(1) PrxdBm(end) -8 0])
%         xlabel('Received power (dBm)', 'FontSize', 12)
% %         xlabel('APD gain (dB)','FontSize',12)
%     end
%     ylabel('log_{10}(BER)', 'FontSize', 12)
%     set(gca, 'FontSize', 12)
% end

%% save results
if sim.save
        % delete large variables
        sim = rmfield(sim, 'f');
        sim = rmfield(sim, 't');
        Tx.Mod = rmfield(Tx.Mod, 'H');    
    try
        save(filename, 'sim','mpam','Tx','Fibers','Rx','ber','PrxdBm');
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