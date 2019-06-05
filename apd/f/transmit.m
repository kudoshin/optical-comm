function [Erx, yd, mpam, Rx] = transmit(dataTX, sim, Tx, Fiber, Apd, Rx, mpam, apdsetting)

% system frequency responses
[HrxPshape,H] = apd_system_received_pulse_shape(mpam, Tx, Fiber, Apd, Rx, sim);

mpam = mpam.norm_levels;
xk0 = mpam.signal(dataTX); % Modulated PAM signal

%% ============================ Preequalization ===============================

if isfield(sim, 'mzm_predistortion') && strcmpi(sim.mzm_predistortion, 'levels') %% Ordinary PAM with predistorted levels
    assert(strcmp(Tx.Mod.type,'MZM'),'predistortion is set but modulator is not MZM')

    % Ajust levels to desired transmitted power and extinction ratio
    mpamPredist = mpam.mzm_predistortion(Tx.Mod.Vswing, Tx.Mod.Vbias, sim.shouldPlot('PAM levels MZM predistortion'));
    xr = mpamPredist.signal(dataTX); % Modulated PAM signal
    xi = 0;
else
    xk = xk0;
    Nf = sim.Mct/Tx.DAC.ros;
    L = length(sim.f);
    DSmask = false(1,L/Nf);
    DSmask(L/2+1-L/(2*Nf):L/2+L/(2*Nf)) = true;
    f = sim.f(DSmask);
    
%     %% compensate APD receiver BW
%     if isfield(sim, 'preemphasis') && sim.preemphasis
%         % frequency response after square-law detection
%         Hrx = H.apd.*H.w.*H.rx.*Rx.eq.Hff(sim.f/mpam.Rs);
%         h = fftshift(real(ifft(ifftshift(Hrx))));
%         h = h(cumsum(abs(h).^2)/sum(abs(h).^2) > 0.001 & cumsum(abs(h).^2)/sum(abs(h).^2) < 0.999);
% 
%         xk = lsqlin(H,xk,-eye(length(xk)),0);
%     end
    
    %% direct detection
    assert(all(xk>=0))
    xk = sqrt(xk);
    xr = real(xk);
    xi = 0;

    %% CD pre-equalization
    if Fiber.L~=0 && isfield(sim, 'CDeq') && sim.CDeq
        Hcdeq = 1./Fiber.Hdisp(f,Tx.Laser.wavelength);
        xk = ifft(fft(xk).*ifftshift(Hcdeq));
        [xr,xi] = clipsc(xk,0);
    end

    %% Predistortion to compensate for MZM non-linear response
    % This predistorts the analog waveform. If sim.mzm_predistortion ==
    % 'levels', then only the levels are predistorted, which is more realistic
    if isfield(sim, 'mzm_predistortion') && strcmpi(sim.mzm_predistortion, 'analog')    
        mzmPredist =@(x) 2/pi*asin(x);
        xr = mzmPredist(xr);
        xi = mzmPredist(xi); % apply predistortion
    end

    xk = xr + 1j*xi;
    %% preemphasis to compensate transmitter and modulator BW
    if isfield(sim, 'preemphasis') && sim.preemphasis
        preemphasis_filter = 1./(H.dac(DSmask).*H.mod(DSmask));
        preemphasis_filter(abs(f)>sim.preemphRange) = 1;

        xk = ifft(fft(xk).*ifftshift(preemphasis_filter));
        xr = real(xk);
        xi = imag(xk);
    end
end
%% DAC
if any(xi~=0)
    xti = dac(xi, Tx.DAC, sim);
else
    xti = xi;
end
xt = dac(xr, Tx.DAC, sim, sim.shouldPlot('DAC output'));
idealdac = Tx.DAC;
idealdac.resolution = Inf;
xt0 = dac(xk0,idealdac,sim);

%% Generate optical signal
Tx.Laser.PdBm = Watt2dBm(Tx.Ptx);
Tx.Laser.H = @(f) Tx.Mod.filt.H(f/sim.fs);
Ecw = Tx.Laser.cw(sim);
if strcmp(Tx.Mod.type, 'MZM')
    Etx = mzm(Ecw, xt, Tx.Mod.filt.H(sim.f/sim.fs))+...
        1j*mzm(Ecw, xti, Tx.Mod.filt.H(sim.f/sim.fs)); % transmitted electric field
    Etx = sqrt(2)*Etx;
elseif strcmp(Tx.Mod.type, 'DML')
    xt = xt - min(xt);
    Etx = Tx.Laser.modulate(xt, sim);
elseif strcmp(Tx.Mod.type, 'EAM')
    if ~isa(Tx.Mod.H,'function_handle')
        Tx.Mod.H =@(f) Tx.Mod.filt.H(f/sim.fs);
    end
    [Etx, ~] = eam(Ecw, 2*xt, Tx.Mod, sim.f);
else
    error('ber_apd_montecarlo: Invalid modulator type. Expecting Tx.Mod.type to be either MZM, EAM, or DML')
end
power = mean(abs(Etx.^2));
%% Ensures that transmitted power is at the right level
Padj = Tx.Ptx/power;
Etx = Etx*sqrt(Padj);

%% Fiber propagation
Erx = Fiber.linear_propagation(Etx, sim.f, Tx.Laser.wavelength);

%% Direct detect
if strcmp(apdsetting,'no noise')
    yt = Apd.detect(Erx, sim.fs, 'no noise');
else
    yt = Apd.detect(Erx, sim.fs, 'gaussian', Rx.N0);
end

%% Whitening filter
if sim.WhiteningFilter
    [~, yt] = Apd.Hwhitening(sim.f, mean(abs(Erx).^2), Rx.N0, yt);
end

%% Automatic gain control
% Normalize signal so that highest level is equal to 1
mpam = mpam.norm_levels;
Rx.AGC = 1/(2*Padj*Tx.Ptx*Apd.Geff*Fiber.link_attenuation(Tx.Laser.wavelength));
yt = yt*Rx.AGC;
yt = yt - mean(yt) + mean(mpam.a);

%% ADC
% ADC performs filtering, quantization, and downsampling
% For an ideal ADC, ADC.ENOB = Inf
% Align received and transmitted signalsclose all
% Rx.ADC.offset = 0;
Rx.ADC.timeRefSignal = xt0; % align filtered signal ytf to this reference
switch lower(Rx.filtering)
    case 'antialiasing' % receiver filter is specified in ADC.filt
        yk = adc(yt, Rx.ADC, sim);
    case 'matched' % receiver filter is matched filter
        Hrx = conj(HrxPshape);
        yk = adc(yt, Rx.ADC, sim, Hrx); % replace ADC antialiasing filter by matched filter
    otherwise
        error('ber_preamp_sys_montecarlo: Rx.filtering must be either antialiasing or matched')
end     

%% Equalization
Rx.eq.trainSeq = dataTX;
[yd, Rx.eq] = equalize(Rx.eq, yk, HrxPshape, mpam, sim);
