function [berenum, pe, distr] = ber_apd_enumeration(mpam, Tx, Fiber, Apd, Rx, sim)
%% Calculate BER for unamplified IM-DD link with APD using enumeration method

%% Pre calculations
Nsymb = mpam.const_size^sim.L; % number of data symbols 
N = sim.Mct*(Nsymb*mpam.symdim + 2*sim.Ndiscard); % total number of points 
    
% change receiver filtering and equalization type
if ~strcmpi(Rx.eq.type, 'none') % if equalization is necessary
    Rx.filtering = 'matched';
    Rx.eq.type = 'Fixed TD-SR-LE'; % always use fixed time-domain symbol rate LE for analysis
    Rx.eq.ros = 1;
end

% Frequency and time

[sim.f, sim.t] = freq_time(N, sim.fs);

% system frequency responses
[HrxPshape, H] = apd_system_received_pulse_shape(mpam, Tx, Fiber, Apd, Rx, sim);


%% Modulated PAM signal

dataTX = debruijn_sequence(mpam.const_size, sim.L).'; % de Bruijin sequence

dataTXext = wextend('1D', 'ppd', dataTX, sim.Ndiscard/mpam.symdim); % periodic extension
dataTXext([1:2*sim.L/mpam.symdim end-2*sim.L/mpam.symdim:end]) = 0; % set 2L first and last symbols to zero


if ~isfield(Tx.Mod,'type')
    Tx.Mod.type = 'EAM';
end


if isfield(sim, 'mzm_predistortion') && strcmpi(sim.mzm_predistortion, 'levels') %% Ordinary PAM with predistorted levels
    assert(strcmp(Tx.Mod.type,'MZM'),'predistortion is set but modulator is not MZM')
    mpam = mpam.norm_levels;
    mpamPredist = mpam.mzm_predistortion(Tx.Mod.Vswing, Tx.Mod.Vbias, sim.shouldPlot('PAM levels MZM predistortion'));
    xk = mpamPredist.signal(dataTXext); % Modulated PAM signal
else
    % Ajust levels to desired extinction ratio
    mpam = mpam.adjust_levels(1,Tx.Mod.rexdB); 
    mpam = mpam.norm_levels;
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
Tx.Laser.PdBm = Watt2dBm(Tx.Ptx);
Ecw = Tx.Laser.cw(sim);
sim.RIN = RIN;

% Modulate
if strcmp(Tx.Mod.type, 'MZM')
    Etx = mzm(Ecw, xt, Tx.Mod.filt.H(sim.f/sim.fs)); % transmitted electric field
    Etx = sqrt(2)*Etx;
elseif strcmp(Tx.Mod.type, 'DML')
    xt = xt - min(xt);
    Etx = Tx.Laser.modulate(xt, sim);
else
    if ~isa(Tx.Mod.H,'function_handle')
        Tx.Mod.H =@(f) Tx.Mod.filt.H(f/sim.fs);
    end
    [Etx, ~] = eam(Ecw, 2*xt, Tx.Mod, sim.f);
end

% Ensures that transmitted power is at the right level
Padj = Tx.Ptx/mean(abs(Etx).^2);
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
AGC = 1/(2*Padj*Tx.Ptx*Rx.PD.Geff*Fiber.link_attenuation(Tx.Laser.wavelength));
yt = yt*AGC;
yt = yt - mean(yt) + mean(mpam.a);

%% ADC
% ADC performs filtering, quantization, and downsampling
% For an ideal ADC, ADC.ENOB = Inf 
Rx.ADC.timeRefSignal = xt; % align filtered signal ytf to this reference
if ~strcmpi(Rx.eq.type, 'none')
    [yk, ~, ~] = adc(yt, Rx.ADC, sim, conj(HrxPshape));
else
    [yk, ~, ~] = adc(yt, Rx.ADC, sim);        
end

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
% if not(mpam.optimize_level_spacing) && isfield(sim, 'mpamOpt') && not(isempty(sim.mpamOpt)) % use threshlds swept in montecarlo simulation
%     Pthresh = zeros(mpam.M-1, 1);
%     mpamOpt = sim.mpamOpt;
%     for k = 1:mpam.M-1
%         Pthresh(k) = (mpam.a(k+1)-mpam.a(k))/(mpamOpt.a(k+1)-mpamOpt.a(k))*(mpamOpt.b(k) - mpamOpt.a(k)) + mpam.a(k);
%     end
% else
%     Pthresh = mpam.b; % decision thresholds referred to the receiver
% end

pe = zeros(mpam.const_size, 1); % symbol error probability for each level

dat = mpam.map(dataTX);

% sort by symbol
indsort = dat == 0:mpam.M-1;
[ysort,varsort] = deal(zeros(max(sum(indsort)),mpam.M));
for m=1:mpam.M
    len = sum(indsort(:,m));
    ysort(1:len,m) = yd(indsort(:,m));
    ysort(len+1:end,m) = nan;
    varsort(1:len,m) = Sshd(indsort(:,m))';
    varsort(len+1:end,m) = nan;
end

distr.condpr =@(x,m) nanmean(normpdf(reshape(x,1,[]),ysort(:,m),sqrt(varsort(:,m))));
distr.tailpr =@(x,m) nanmean(qfunc((reshape(x,1,[])-ysort(:,m))./sqrt(varsort(:,m))));
distr.mean = nanmean(ysort);
distr.var = nanmean(varsort)+nanvar(ysort);

% use tailpr to optimize thresholds
mpam = mpam.optimize_thresholds_enum(distr);

if mpam.M ~=3
    pe(1) = distr.tailpr(mpam.b(1),1);
    pe(2:mpam.M-1) = distr.tailpr(mpam.b(2:mpam.M-1)',2:mpam.M-1) + ...
        1 - distr.tailpr(mpam.b(1:mpam.M-2)',2:mpam.M-1);
    pe(mpam.M) = 1 - distr.tailpr(mpam.b(mpam.M-1),mpam.M);
else
    bit_assign = [0 2 6 1 3 7 4 5];
    d = 3*squareform(pdist(de2bi(bit_assign),'hamming'));
    
    Qr = distr.tailpr(mpam.b,1:mpam.M-1);   % right xover pr.
    Ql = 1- distr.tailpr(mpam.b,2:mpam.M); % left xover pr.
    
    x = linspace(mpam.b(2),3*mpam.a(end),500)';
    p68 = Ql(2)*Qr(2) + trapz(x,distr.condpr(x,3).*distr.tailpr(x,2));
    
    pe(1) = (Qr(1)-Qr(1)^2)*(d(1,2)+d(1,4))+Qr(1)^2*d(1,5);
    pe(2) = (Ql(1)-Ql(1)*Qr(1))*d(1,2)+(Qr(2)-Qr(2)*Qr(1))*d(2,3)+...
        (Qr(1)-Qr(1)*Ql(1)-Qr(1)*Qr(2))*d(2,5)+Ql(1)*Qr(1)*d(2,4)+...
        Qr(2)*Qr(1)*d(2,6);
    pe(3) = (Ql(2)-Ql(2)*Qr(1))*d(2,3) + Ql(2)*Qr(1)*d(3,5) +...
        (Qr(1)-Qr(1)*Ql(2))*d(3,6);
    pe(4) = (Ql(1)-Ql(1)*Qr(1))*d(1,4)+Qr(1)*Ql(1)*d(2,4)+...
        (Qr(1)-Qr(1)*Ql(1)-Qr(1)*Qr(2))*d(4,5)+...
        (Qr(2)-Qr(2)*Qr(1))*d(4,7)+Qr(1)*Qr(2)*d(2,8);
    pe(5) = Ql(1)^2*d(1,5)+Ql(1)*(1-Ql(1)-Qr(2))*d(2,5)+...
        Ql(1)*Qr(2)*d(3,5)+Ql(1)*(1-Ql(1)-Qr(2))*d(4,5)+...
        Qr(2)*(1-Ql(1)-Qr(2)/2)*d(5,6)+Ql(1)*Qr(2)*d(5,7)+...
        Qr(2)*(1-Ql(1)-Qr(2)/2)*d(5,8);
    pe(6) = Ql(2)*Ql(1)*d(2,6)+Ql(1)*(1-Ql(2))*d(3,6)+...
        Ql(2)*(1-Ql(1)-Qr(2))*d(5,6)+p68*d(6,8);
    pe(7) = Ql(2)*(1-Qr(1))*d(4,7)+Qr(1)*Ql(2)*d(5,7)+...
        Qr(1)*(1-Ql(2))*d(7,8);
    pe(8) = Ql(1)*Ql(2)*d(4,8)+Ql(2)*(1-Ql(1)-Qr(2))*d(5,8)+p68*d(6,8)+...
        Ql(1)*(1-Ql(2))*d(7,8);
end
berenum = mean(pe)/log2(mpam.const_size); 

if sim.shouldPlot('Conditional PDF') 
    figure(110)
    set(gca,'colororderindex',1)
    for m=1:mpam.M
        x=linspace(distr.mean(m)-4*sqrt(distr.var(m)),...
            distr.mean(m)+4*sqrt(distr.var(m)));
        plot(x,distr.condpr(x,m));
        hold on;
    end
    grid on
    xlabel('normalized voltage')
    ylabel('pdf')
%     set(gca,'colororderindex',1)
%     plot(linspace(0,1),normpdf(linspace(0,1)',distr.mean,sqrt(distr.var)),'--')
    vline(mpam.b)
    hold off
    title(num2str(berenum,'BER = %.2e'))
    drawnow
end
