function [ber, dataTX, Etx] = ber_apd_montecarlo(mpam, Tx, Fiber, Apd, Rx, sim, dataTX, Etx)
%% Calculate BER of unamplified IM-DD system with APD detector using montecarlo simulation

% Ajust levels to desired transmitted power and extinction ratio
% mpam = mpam.adjust_levels(Tx.Ptx, Tx.Mod.rexdB);

if ~exist('dataTX','var') || isempty(dataTX)
    %% Modulated PAM signal
    % dataTX = [0 0 0 0 0 1 0 0 0 0 0];
    dataTX= randi([0 mpam.const_size-1], 1, sim.Nsymb/mpam.symdim); % Random sequence
    Nzero = 10;
    dataTX([1:Nzero/mpam.symdim end-Nzero/mpam.symdim+1:end]) = 0;

    Etx = transmit(dataTX, sim, Tx, Fiber, Apd, Rx, mpam);
end
[~, yd, Rx] = receive(dataTX, sim, Tx, Fiber, Apd, Rx, mpam, Etx,'gaussian');  

% Symbols to be discard in BER calculation
dataTXcut = dataTX;
ndiscard = [1:Rx.eq.Ndiscard(1)+sim.Ndiscard (sim.Nsymb-Rx.eq.Ndiscard(2)-sim.Ndiscard+1):sim.Nsymb];
ydfull = yd;
yd(ndiscard) = []; 
dataTXcut(round(ndiscard(mpam.symdim:mpam.symdim:end)/mpam.symdim)) = [];

%% Demodulate
if mpam.optimize_level_spacing
    dataRX = mpam.demod(yd);
else
    [dataRX, mpam] = mpam.demod_sweeping_thresholds(yd, dataTXcut);
end

%% True BER
[~, ber] = biterr(dataRX, dataTXcut);

%% Plots
if sim.shouldPlot('Empirical noise pdf')
    % Empirical pdf for a level
    figure(100)
    [nn, xx] = hist(yd(dataTXcut == 2), 50);
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
    isort = mpam.map(dataTXcut)==(1:mpam.M)'-1;
    for m=1:mpam.M
        [nn, xx] = hist(yd(isort(m,:)), 50);
        nn = nn/trapz(xx, nn);
        plot(xx, nn,'-.')
        hold on
    end
    hold off
%     vline(mpam.b,'b:')
    drawnow
end
