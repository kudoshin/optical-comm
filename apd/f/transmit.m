function [Erx, yd, mpam, Rx] = transmit(dataTX, sim, Tx, Fiber, Apd, Rx, mpam, apdsetting, silent)

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
    L = length(sim.f);
    Nf = sim.Mct/Tx.DAC.ros;
    DSmask = false(1,L);
    DSmask(L/2+1-L/(2*Nf):L/2+L/(2*Nf)) = true;
    f = sim.f(DSmask);
    
    %% compensate APD receiver BW
    if isfield(sim, 'preemphasis') && sim.preemphasis
        % frequency response from APD after square-law detection
        Hrx = H.apd(DSmask).*H.w(DSmask);
        
        Heq = 1./Hrx;
        Heq(abs(f)>mpam.Rs) = 1;
        
        xk = real(ifft(ifftshift(Heq).*fft(xk0)));
        xk(xk<0) = 0;
        
%         % trying to do a minimization
%         h = fftshift(real(ifft(ifftshift(Hrx))));
%         h = h(cumsum(abs(h).^2)/sum(abs(h).^2) > 0.001 & cumsum(abs(h).^2)/sum(abs(h).^2) < 0.999);
% 
%         h0 = round(grpdelay(h, 1, 1))+1;
%         hd_prev = h(h0:-Nf:1);
%         hd_post = h(h0:Nf:end);
%         hd = [hd_prev(end:-1:2) hd_post];
%         hd_prev = hd_prev/sum(abs(hd).^2);
%         hd_post = hd_post/sum(abs(hd).^2);
%         
%         H = toeplitz([hd_post zeros(1,L/Nf-length(hd_post))],...
%             [hd_prev zeros(1,L/Nf-length(hd_prev))]);
%         xk = lsqnonneg(H,xk);

    end
    
    %% direct detection
    assert(all(xk>=0))
    xk = sqrt(xk);
    xr = real(xk);
    xi = 0;

    %% CD pre-equalization
    if isfield(sim, 'CDeq') && sim.CDeq
        CDFiber = Fiber;
        CDFiber.L = 25e3;
        Hcdeq = 1./CDFiber.Hdisp(f,Tx.Laser.wavelength);
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
    %% preemphasis to compensate modulator BW
    if isfield(sim, 'preemphasis') && sim.preemphasis
        preemphasis_filter = 1./H.mod(DSmask);
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
    xti = 0;
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
        [yk,~,yt] = adc(yt, Rx.ADC, sim);
    case 'matched' % receiver filter is matched filter
        Hrx = conj(HrxPshape);
        [yk,~,yt] = adc(yt, Rx.ADC, sim, Hrx); % replace ADC antialiasing filter by matched filter
    otherwise
        error('ber_preamp_sys_montecarlo: Rx.filtering must be either antialiasing or matched')
end     

%% Equalization
Rx.eq.trainSeq = dataTX;
[yd, Rx.eq] = equalize(Rx.eq, yk, HrxPshape, mpam, sim);

if ~exist('silent','var') || ~silent
    if sim.shouldPlot('Optical eye diagram')  
        Ntraces = 500;
        Nstart = sim.Ndiscard*sim.Mct + 1;
        Nend = min(Nstart + Ntraces*2*sim.Mct, length(Etx));
        figure(103), clf, box on
        eyediagram(abs(Etx(Nstart:Nend)).^2, 2*sim.Mct)
        title('Optical eye diagram')
        drawnow
    end

    if sim.shouldPlot('Received signal eye diagram')
        Ntraces = 500;
        Nstart = sim.Ndiscard*sim.Mct + 1;
        Nend = min(Nstart + Ntraces*2*sim.Mct, length(Erx));
        figure(104), clf, box on
        eyediagram(yt(Nstart:Nend), 2*sim.Mct)
        title('Noiseless received signal eye diagram')
        hold on
        h1 = hline(mpam.a, '-k');
        h2 = hline(mpam.b, '--k');
        h3 = vline(sim.Mct+1, 'k');
        legend([h1(1) h2(1) h3], {'Levels', 'Decision thresholds', 'Sampling point'})
        drawnow
    end
end

