function [ber_awgn, noise_std]  = ber_apd_awgn(mpam, Tx, Fiber, Apd, Rx, sim)
%% Calculate BER for unamplified IM-DD link with APD using AWGN Approximation including noise enhancement penalty
% change receiver filtering and equalization type
if ~strcmpi(Rx.eq.type, 'none') % if equalization is necessary
    Rx.filtering = 'matched';
    Rx.eq.type = 'Fixed TD-SR-LE'; % always use fixed time-domain symbol rate LE for analysis
    Rx.eq.ros = 1;
end

% system frequency responses
[HrxPshape, H] = apd_system_received_pulse_shape(mpam, Tx, Fiber, Apd, Rx, sim);

% Equalization
[~, Rx.eq] = equalize(Rx.eq, [], HrxPshape, mpam, sim);

% Noise std: includes RIN, shot and thermal noise (assumes gaussian stats)
noiseSTD = Apd.stdNoise(H.w.*H.rx, Rx.eq.Hff(sim.f/(Rx.eq.ros*mpam.Rs)), Rx.N0, Tx.Laser.RIN, sim);

% Receiver filter
% For symbol-rate sampling linear equalizer = APD -> Whitening filter ->
% matched filter -> equalizer (in continuous time)
Hshot = H.apd.*H.w.*H.rx.*Rx.eq.Hff(sim.f/mpam.Rs); % Shot noise shape

h2 = fftshift(real(ifft(ifftshift(Hshot))));
h = h2(cumsum(abs(h2).^2)/sum(abs(h2).^2) > 0.001 & cumsum(abs(h2).^2)/sum(abs(h2).^2) < 0.999);
hh = h.*conj(h); % |h(t)|^2
hh = hh/abs(sum(hh)); % normalize

hh0 = round(grpdelay(hh, 1, 1));
hhd_prev = hh(hh0:-sim.Mct:1);
hhd_post = hh(hh0:sim.Mct:end);
hhd = [hhd_prev(end:-1:2) hhd_post];
hh0 = length(hhd_prev);
hhd = hhd/abs(sum(hhd));
isi_sh = sum(hhd) - hhd(hh0);

% if strcmpi(Rx.eq.type, 'none')
%     % uncorrected ISI
%     Hch = H.pshape.*H.ch.*H.adc;
%     h1 = fftshift(real(ifft(ifftshift(Hch))));
%     hh1 = h1.*conj(h1);
%     H2 = fft(ifftshift(hh1)).*ifftshift(H.adc);
%     h2 = fftshift(real(ifft(H2)));
%     hh = h2(cumsum(abs(h1).^2)/sum(abs(h1).^2) > 0.001 & cumsum(abs(h1).^2)/sum(abs(h1).^2) < 0.999);
% %     hh = h.*conj(h); % |h(t)|^2
% 
%     hh0 = round(grpdelay(hh, 1, 1));
%     hhd_prev = hh(hh0:-sim.Mct:1);
%     hhd_post = hh(hh0:sim.Mct:end);
%     hhd = [hhd_prev(end:-1:2) hhd_post];
%     hh0 = length(hhd_prev);
%     hhd = hhd/sum(hhd);
%     isi = sum(hhd)-hhd(hh0);
% else
%     hhd = 1;
%     hh0 = 1;
%     isi = 0;
% end

% Refer PAM levels to the receiver
mpam = mpam.adjust_levels(Apd.Geff*Tx.Ptx*Fiber.link_attenuation(Tx.Laser.wavelength), Tx.Mod.rexdB);

% Psymbs = mean(mpam.a); % symbol power in the memory of APD
Psymbs = mpam.a(end);

% Update noise standard deviation at each level including Shot noise ISI
% noise_std = @(P) sqrt(noiseSTD((1-isi_sh)*P*hhd(hh0).^2  + Psymbs*isi_sh).^2 + var(mpam.a)*isi);
noise_std = @(P) noiseSTD((1-isi_sh)*P+ Psymbs*isi_sh);
% In calculating noise_std for AWGN channel, all other symbols in hh
% are assumed to be the highest PAM level.

% if not(mpam.optimize_level_spacing)
    mpam = mpam.optimize_thresholds(noise_std);
% end

% Calculate BER
ber_awgn = mpam.berAWGN(@(P) noise_std(P));
