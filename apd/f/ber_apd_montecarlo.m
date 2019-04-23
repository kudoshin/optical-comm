function [ber, mpam] = ber_apd_montecarlo(mpam, Tx, Fiber, Apd, Rx, sim)
%% Calculate BER of unamplified IM-DD system with APD detector using montecarlo simulation

% System received pulse shape frequency response
HrxPshape = apd_system_received_pulse_shape(mpam, Tx, Fiber, Apd, Rx, sim); % this is only used if rx.filtering = matched or eq.type = fixed...

% Ajust levels to desired transmitted power and extinction ratio
% mpam = mpam.adjust_levels(Tx.Ptx, Tx.Mod.rexdB);

%% Modulated PAM signal
% dataTX = [0 0 0 0 0 1 0 0 0 0 0];
if mpam.M ~=3
    dataTX= randi([0 mpam.M-1], 1, sim.Nsymb); % Random sequence
else
    dataTX = randi([0 7], 1, ceil(sim.Nsymb/2));
end
Nzero = 10;
dataTX([1:Nzero/mpam.symdim end-Nzero/mpam.symdim+1:end]) = 0;

if isfield(sim, 'precomp') && sim.precomp
    mzmtf =@(v) sin(pi/2*v);
    mzmitf =@(v) 2/pi*asin(v);
    pmax = mzmtf(Tx.Mod.Vswing/2+Tx.Mod.Vbias)^2;
    pmin = mzmtf(-Tx.Mod.Vswing/2+Tx.Mod.Vbias)^2;
    mpampre = mpam.norm_levels();
    mpampre = mpampre.set_levels(mpampre.a*(pmax-pmin)+pmin,...
        mpampre.b*(pmax-pmin)+pmin);
    x0 = mpampre.signal(dataTX);   % desired PAM signal
%     % preemphasis for APD/ADC filtering
%     Hpre = 1./Rx.ADC.filt.H(sim.f/sim.fs);
%     Hpre(abs(sim.f)>18e9)=0;
%     x1 = ifft(fft(x0).*ifftshift(Hpre)); 
    x1=x0;
    x2 = sqrt(x1);
    % undo dispersion
    Fiber.L = -Fiber.L;
    x3 = Fiber.linear_propagation(x2,sim.f,Tx.Laser.wavelength);
    Fiber.L = -Fiber.L;
    % clip
    x4=clip(x3,[0 1]);
    % undo MZM nonlinearity
    x5 = mzmitf(real(x4))+1j*mzmitf(imag(x4));
    % preemphasis for DAC/MZM filtering
    Hpre = 1./(Tx.DAC.filt.H(sim.f/sim.fs).*Tx.Mod.filt.H(sim.f/sim.fs));
    Hpre(abs(sim.f)>50e9)=0;
    xk = ifft(fft(x5).*ifftshift(Hpre));
%     xk=x5;
    xt = dac(real(xk),Tx.DAC,sim, sim.shouldPlot('DAC output')) + 1j*dac(imag(xk),Tx.DAC,sim);
else
    if isfield(sim, 'mzm_predistortion') && strcmpi(sim.mzm_predistortion, 'levels') %% Ordinary PAM with predistorted levels
        assert(strcmp(Tx.Mod.type,'MZM'),'predistortion is set but modulator is not MZM')

        % Ajust levels to desired transmitted power and extinction ratio
        mpam=mpam.norm_levels;
        mpamPredist = mpam.mzm_predistortion(Tx.Mod.Vswing, Tx.Mod.Vbias, sim.shouldPlot('PAM levels MZM predistortion'));
        xk = mpamPredist.signal(dataTX); % Modulated PAM signal
    else
        mpam = mpam.adjust_levels(1,Tx.Mod.rexdB);
        mpam = mpam.norm_levels;
        xk = mpam.signal(dataTX); % Modulated PAM signal
    end  

    %% DAC
    xt = dac(xk, Tx.DAC, sim, sim.shouldPlot('DAC output')); 
end

%% Generate optical signal
Tx.Laser.PdBm = Watt2dBm(Tx.Ptx);
Tx.Laser.H = @(f) Tx.Mod.filt.H(f/sim.fs);
Ecw = Tx.Laser.cw(sim);
if strcmp(Tx.Mod.type, 'MZM')
    Etx = mzm(Ecw, xt, Tx.Mod.filt.H(sim.f/sim.fs)); % transmitted electric field
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

%% Detect and add noises
% yt = Apd.detect(Erx, sim.fs, 'no noise');
yt = Apd.detect(Erx, sim.fs, 'gaussian', Rx.N0);

%% Whitening filter
if sim.WhiteningFilter
    [~, yt] = Apd.Hwhitening(sim.f, mean(abs(Erx).^2), Rx.N0, yt);
end

%% Automatic gain control
% Normalize signal so that highest level is equal to 1
mpam = mpam.norm_levels;
AGC = 1/(2*Padj*Tx.Ptx*Apd.Geff*Fiber.link_attenuation(Tx.Laser.wavelength));
yt = yt*AGC;
yt = yt - mean(yt) + mean(mpam.a);

%% ADC
% ADC performs filtering, quantization, and downsampling
% For an ideal ADC, ADC.ENOB = Inf
% Align received and transmitted signals
% Rx.ADC.offset = 0;
Rx.ADC.timeRefSignal = xt; % align filtered signal ytf to this reference
switch lower(Rx.filtering)
    case 'antialiasing' % receiver filter is specified in ADC.filt
        [yk, ~, ytf] = adc(yt, Rx.ADC, sim);
    case 'matched' % receiver filter is matched filter
        Hrx = conj(HrxPshape);
        [yk, ~, ytf] = adc(yt, Rx.ADC, sim, Hrx); % replace ADC antialiasing filter by matched filter
    otherwise
        error('ber_preamp_sys_montecarlo: Rx.filtering must be either antialiasing or matched')
end     

%% Equalization
Rx.eq.trainSeq = dataTX;
[yd, Rx.eq] = equalize(Rx.eq, yk, HrxPshape, mpam, sim);

% Symbols to be discard in BER calculation
ndiscard = [1:Rx.eq.Ndiscard(1)+sim.Ndiscard (sim.Nsymb-Rx.eq.Ndiscard(2)-sim.Ndiscard+1):sim.Nsymb];
ydfull = yd;
yd(ndiscard) = []; 
dataTX(round(ndiscard(mpam.symdim:mpam.symdim:end)/mpam.symdim)) = [];

%% Demodulate
% if mpam.optimize_level_spacing
%     dataRX = mpam.demod(yd);
% else
    [dataRX, mpam] = mpam.demod_sweeping_thresholds(yd, dataTX);
% end

%% True BER
[~, ber] = biterr(dataRX, dataTX);

%% Plots
if sim.shouldPlot('Empirical noise pdf')
    % Empirical pdf for a level
    figure(100)
    [nn, xx] = hist(yd(dataTX == 2), 50);
    nn = nn/trapz(xx, nn);
    bar(xx, nn)
    title('Empirical pdf for PAM level 2')
end    

if sim.shouldPlot('Equalizer')
    figure(101), clf
    for k = 1:size(Rx.eq.h, 2)
        [h, w] = freqz(Rx.eq.h(:, k), 1);
        subplot(121), hold on, box on
        plot(w/(2*pi), abs(h).^2)
        
        subplot(122), hold on, box on
        plot(w/(2*pi), unwrap(angle(h)))
    end
    subplot(121), hold on, box on
    xlabel('Normalized frequency')
    ylabel('|H(f)|^2')
    title('Equalizer amplitude response')

    subplot(122), hold on, box on
    xlabel('Normalized frequency')
    ylabel('arg(H(f))')
    title('Equalizer phase response')
    drawnow
end

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
    eyediagram(abs(Erx(Nstart:Nend)).^2, 2*sim.Mct)
    title('Received optical eye diagram (before noise)')
    mpam = mpam.norm_levels();
    Ntraces = 500;
    Nstart = sim.Ndiscard*sim.Mct + 1;
    Nend = min(Nstart + Ntraces*2*sim.Mct, length(ytf));
    figure(105), clf, box on, hold on
    eyediagram(ytf(Nstart:Nend), 2*sim.Mct)
    title('Received signal eye diagram')
    a = axis;
    h1 = plot(a(1:2), mpam.a*[1 1], '-k');
    h2 = plot(a(1:2), mpam.b*[1 1], '--k');
    h3 = plot((sim.Mct+1)*[1 1], a(3:4), 'k');
    legend([h1(1) h2(1) h3], {'Levels', 'Decision thresholds', 'Sampling point'})
    drawnow
end

if sim.shouldPlot('Signal after equalization')
    mpam = mpam.norm_levels();
    figure(106), clf, box on, hold on
    h1 = plot(ydfull, 'o');
    a = axis;
    h2= plot(a(1:2), mpam.a*[1 1], '-k');
    h3 = plot(a(1:2), mpam.b*[1 1], '--k');
    h4 = plot(Rx.eq.Ndiscard(1)*[1 1], a(3:4), ':k');
    h5 = plot((sim.Nsymb-Rx.eq.Ndiscard(2))*[1 1], a(3:4), ':k');
    legend([h1 h2(1) h3(1) h4], {'Equalized samples', 'PAM levels',...
        'Decision thresholds', 'BER measurement window'})
    title('Signal after equalization')
    axis([1 sim.Nsymb -0.2 1.2])
    drawnow
end

if sim.shouldPlot('Conditional PDF')
    figure(110)
    set(gca,'colororderindex',1)
    isort = mpam.map(dataTX)==(1:mpam.M)'-1;
    for m=1:mpam.M
        [nn, xx] = hist(yd(isort(m,:)), 50);
        nn = nn/trapz(xx, nn);
        plot(xx, nn,'-.')
        hold on
    end
    hold off
    drawnow
end
