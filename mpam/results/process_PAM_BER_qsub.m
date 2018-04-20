%% Process data saved by PAM_BER_qsub.m
clear, clc, close all

addpath ../
addpath ../../f/
addpath ../../apd/
addpath ../../apd/f/

M = [2 4 8];
folder = '../../access_nets/results/4PAM/';

% PAM_BER_L=0.5km_lamb=1380nm_ModBW=30GHz_amplified=0_Ntaps=9_ENOB=5_ros=1.25

BERtarget = 1.8e-4;
Amplified = {'0', '1', 'apd'};
amplified = {'pin','soa','apd'};
ModBWGHz = 10;
ENOB = 5;
Ntaps = 15;
lamb = 1550;
ros = 2;

LineStyle = {'-', '--'};
Marker = {'o', 's', 'v'};
Color = {[51, 105, 232]/255, [153,153,155]/255, [255,127,0]/255};
% Lkm = 0:0.5:10;
Lkm = [2 10 20 40 50 55];

% Fiber = fiber();
D = zeros(1, length(Lkm));
PrxdBm = zeros(length(Amplified), length(Lkm));
% for a = 1:length(Amplified)
%     for k = 1:length(Lkm)
for a = 1:3
    for k=1:length(Lkm)  
        filename = [folder sprintf('L=%skm/PAM_BER_lamb=%dnm_ModBW=%dGHz_amplified=%s_RecBW=Inf_Ntaps=%d_ENOB=%d_ros=%.2f.mat',...
            num2str(Lkm(k)), lamb, ModBWGHz, Amplified{a}, Ntaps, ENOB, ros)];  

        try 
            S = load(filename, '-mat');
            D(k) = 1e6*S.Fibers(1).D(S.Tx.Laser.wavelength)*S.Fibers(1).L/1e3;

            % Realizations were already averaged
            BERcount = log10(S.ber.count);
            BERgauss = log10(S.ber.gauss);

            idx = find(BERcount <= -2 & BERcount >= -5);
            [PrxdBm(a, k), f] = fit_ber(S.Tx.PtxdBm(idx), S.ber.count(idx), BERtarget);
            
            % Plot
            figure(2), clf, hold on, box on
            hline = plot(S.Tx.PtxdBm, BERcount, '-o');
            plot(S.Tx.PtxdBm, f(S.Tx.PtxdBm), '-', 'Color', get(hline, 'Color'));
            plot(S.Tx.PtxdBm, BERgauss, '--', 'Color', get(hline, 'Color'));
            axis([S.Tx.PtxdBm([1 end]) -8 0])
            title(sprintf('L = %.1f km, D = %.2f', S.Fibers(1).L/1e3, D(k)))
            drawnow   
        catch e
            filename
            warning(e.message)
            PrxdBm(a, k) = NaN;
        end
    end
    figure(1), hold on, box on
    plot(D, PrxdBm(a, :), '-o', 'Color', Color{a}, 'LineWidth', 2,...
        'MarkerFaceColor', 'w', 'DisplayName', sprintf(amplified{a}));
end

figure(1)
xlabel('Dispersion (ps/nm)')
ylabel('Receiver sensitivity (dBm)')
legend('-dynamiclegend')
