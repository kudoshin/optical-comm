function [Erx, yd, Rx] = receive(dataTX, sim, Tx, Fiber, Apd, Rx, mpam, Etx, apdsetting, silent)

% system frequency response
HrxPshape = apd_system_received_pulse_shape(mpam, Tx, Fiber, Apd, Rx, sim);

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
idealdac = Tx.DAC;
idealdac.resolution = Inf;
xk0 = mpam.signal(dataTX); % Modulated PAM signal
xt0 = dac(xk0,idealdac,sim);
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

