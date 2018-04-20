function berenum = ber_apd_enumeration(mpam, Tx, Fiber, Apd, Rx, sim)
%% Calculate BER for unamplified IM-DD link with APD using enumeration method

%% Pre calculations
Nsymb = mpam.M^sim.L; % number of data symbols 
N = sim.Mct*(Nsymb + 2*sim.Ndiscard); % total number of points 

% change receiver filtering and equalization type
Rx.filtering = 'matched';
if ~strcmpi(Rx.eq.type, 'none') % if equalization is necessary
    Rx.eq.type = 'Fixed TD-SR-LE'; % always use fixed time-domain symbol rate LE for analysis
    Rx.eq.ros = 1;
end

% Frequency and time

[sim.f, sim.t] = freq_time(N, sim.fs);

% system frequency responses
[HrxPshape, H] = apd_system_received_pulse_shape(mpam, Tx, Fiber, Apd, Rx, sim);


%% Modulated PAM signal
dataTX = debruijn_sequence(mpam.M, sim.L).'; % de Bruijin sequence
dataTXext = wextend('1D', 'ppd', dataTX, sim.Ndiscard); % periodic extension
dataTXext([1:2*sim.L end-2*sim.L:end]) = 0; % set 2L first and last symbols to zero

if ~isfield(Tx.Mod,'type')
    Tx.Mod.type = 'EAM';
end
% predistortion assumes equal level spacing
if isfield(sim, 'mzm_predistortion') && strcmpi(sim.mzm_predistortion, 'levels') %% Ordinary PAM with predistorted levels
    assert(strcmp(Tx.Mod.type,'MZM'),'predistortion is set but modulator is not MZM')
    mpamPredist = mpam.mzm_predistortion(Tx.Mod.Vswing, Tx.Mod.Vbias, sim.shouldPlot('PAM levels MZM predistortion'));
    xk = mpamPredist.signal(dataTXext); % Modulated PAM signal
    AGC = 1/(Apd.Geff*Fiber.link_attenuation(Tx.Laser.wavelength));
else
    % Ajust levels to desired transmitted power and extinction ratio
    mpam = mpam.adjust_levels(Tx.Ptx, Tx.Mod.rexdB);
    AGC = 1/(mpam.a(end)*Apd.Geff*Fiber.link_attenuation(Tx.Laser.wavelength));
    
    xk = mpam.signal(dataTXext); % Modulated PAM signal
end  

%% ============================ Preemphasis ===============================
if isfield(sim, 'preemphasis') && sim.preemphasis
    femph = abs(freq_time(sim.Nsymb*sim.ros.txDSP, mpam.Rs*sim.ros.txDSP));
    femph(femph >= sim.preemphRange) = 0;
    preemphasis_filter = 10.^(polyval([-0.0013 0.5846 0], femph/1e9)/20);  % Coefficients were measured in the lab  

    xk = real(ifft(fft(xk).*ifftshift(preemphasis_filter)));
end

%% Predistortion to compensate for MZM non-linear response
% This predistorts the analog waveform. If sim.mzm_predistortion ==
% 'levels', then only the levels are predistorted, which is more realistic
if isfield(sim, 'mzm_predistortion') && strcmpi(sim.mzm_predistortion, 'analog')    
    xk = 2/pi*asin(sqrt(abs(xk))).*sign(xk); % apply predistortion
end

%% DAC
xt = dac(xk, Tx.DAC, sim); 

%% Driver
% Adjust gain to compensate for preemphasis
xt = Tx.Vgain*(xt - mean(xt)) + Tx.VbiasAdj*mean(xt);

%% Generate optical signal
if ~isfield(sim, 'RIN')
    sim.RIN = false;
end

RIN = sim.RIN;
sim.RIN = false; % RIN is not modeled here since number of samples is not high enough to get accurate statistics
sim.phase_noise = false; 
Ecw = Tx.Laser.cw(sim);
sim.RIN = RIN;

% Modulate
if strcmp(Tx.Mod.type, 'MZM')
    Etx = mzm(Ecw, xt, Tx.Mod.filt.H(sim.f/sim.fs)); % transmitted electric field
elseif strcmp(Tx.Mod.type, 'DML')
    xt = xt - min(xt);
    Etx = Tx.Laser.modulate(xt, sim);
else
    if ~isa(Tx.Mod.H,'function_handle')
        Tx.Mod.H =@(f) Tx.Mod.filt.H(f/sim.fs);
    end
    [Etx, ~] = eam(Ecw, xt, Tx.Mod, sim.f);
end

% Ensures that transmitted power is at the right level
Padj = Tx.Ptx/mean(abs(Etx).^2);
AGC = AGC/Padj;
Etx = Etx*sqrt(Padj);

%% Fiber propagation
Erx = Fiber.linear_propagation(Etx, sim.f, Tx.Laser.wavelength);

%% Direct detect
yt = Apd.detect(Erx, sim.fs, 'no noise');

%% Noise whitening filter
if sim.WhiteningFilter
    [~, yt] = Apd.Hwhitening(sim.f, mean(abs(Erx).^2), Rx.N0, yt);
end

%% Automatic gain control
mpam = mpam.norm_levels;
yt = yt*AGC;
yt = yt - mean(yt) + mean(mpam.a);

%% ADC
% ADC performs filtering, quantization, and downsampling
% For an ideal ADC, ADC.ENOB = Inf
Rx.ADC.ros = 1; % symbol-rate sampling 
Rx.ADC.timeRefSignal = xt; % align filtered signal ytf to this reference
[yk, ~, ytf] = adc(yt, Rx.ADC, sim, conj(HrxPshape));

%% Equalization
[yd, Rx.eq] = equalize(Rx.eq, yk, HrxPshape, mpam, sim);

% Symbols to be discard in BER calculation
yd = yd(sim.Ndiscard+1:end-sim.Ndiscard);

%% Calculate signal-dependent noise variance after matched filtering and equalizer 
Ssh = Apd.varShot(abs(Erx).^2, 1)/2; % two-sided shot noise PSD

% Receiver filter
% For symbol-rate sampling linear equalizer = APD -> Whitening filter ->
% matched filter -> equalizer (in continuous time)
Hshot = H.apd.*H.w.*H.rx.*Rx.eq.Hff(sim.f/mpam.Rs); % Shot noise shape
BWshot = trapz(sim.f, abs(Hshot).^2); % two-sided shot noise bandwidth

Hthermal = H.w.*H.rx.*Rx.eq.Hff(sim.f/mpam.Rs); % thermal noise shape
BWthermal = trapz(sim.f, abs(Hthermal).^2); % two-sided thermal noise bandwidth

h2 = fftshift(real(ifft(ifftshift(Hshot))));
h = h2(cumsum(abs(h2).^2)/sum(abs(h2).^2) > 0.001 & cumsum(abs(h2).^2)/sum(abs(h2).^2) < 0.999);
hh = h.*conj(h); % |h(t)|^2
hh = hh/abs(sum(hh)); % normalize

Sshf = BWshot*conv(hh, Ssh);
Sshf = delay_signal(Sshf, -grpdelay(hh, 1, 1)); % remove delay due to equalizer

% Add thermal noise
Sshf = Sshf + Rx.N0/2*BWthermal;

% Normalize and sample
Sshd = Sshf(1:sim.Mct:end);
Sshd = (AGC)^2*Sshd(sim.Ndiscard+1:end-sim.Ndiscard);

%% Calculate error probabilities using Gaussian approximation for each transmitted symbol
if not(mpam.optimize_level_spacing) && isfield(sim, 'mpamOpt') && not(isempty(sim.mpamOpt)) % use threshlds swept in montecarlo simulation
    Pthresh = zeros(mpam.M-1, 1);
    mpamOpt = sim.mpamOpt;
    for k = 1:mpam.M-1
        Pthresh(k) = (mpam.a(k+1)-mpam.a(k))/(mpamOpt.a(k+1)-mpamOpt.a(k))*(mpamOpt.b(k) - mpamOpt.a(k)) + mpam.a(k);
    end
else
    Pthresh = mpam.b; % decision thresholds referred to the receiver
end

pe = zeros(mpam.M, 1); % symbol error probability for each level
dat = gray2bin(dataTX, 'pam', mpam.M); % fix index mapping
for k = 1:Nsymb
    sig = sqrt(Sshd(k));
    
    if dat(k) == mpam.M-1
        pe(dat(k)+1) = pe(dat(k)+1) + qfunc((yd(k)-Pthresh(end))/sig);
    elseif dat(k) == 0
        pe(dat(k)+1) = pe(dat(k)+1) + qfunc((Pthresh(1)-yd(k))/sig);
    else 
        pe(dat(k)+1) = pe(dat(k)+1) + qfunc((Pthresh(dat(k) + 1) - yd(k))/sig);
        pe(dat(k)+1) = pe(dat(k)+1) + qfunc((yd(k) - Pthresh(dat(k)))/sig);
    end
end

pe = real(pe)/Nsymb;

berenum = sum(pe)/log2(mpam.M);
